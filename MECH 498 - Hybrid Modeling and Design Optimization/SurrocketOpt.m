function [I_sp, max_alt, vel_off_rail, p] = SurrocketOpt(x, trade)
%SURROCKET Models a hybrid rocket from design to apogee
%   'x' is a design vector containing the following variables:
%       'V_tank' - oxidizer tank volume (m^3)
%       'C_inj' - effective injector area (m^2)
%       'L' - fuel grain length (m)
%       'd_port_init' - initial fuel grain port diameter (m)
%       'd_th' - nozzle throat diameter (m)
%       'A_ratio_nozzle' - nozzle area ratio as a fraction
%   'trade' is the index of the parameter to be varied in the trade study
%       i.e. if you wanted to vary fuel grain length, trade = 3

%% INPUTS %%
% Load input files
CEA = load('MuleSim3CEA.mat'); % Combustion product lookup table
N2Osat = load('N2Osat.mat'); % Nitrous oxide saturation properties
CdA = load('Cdavg.mat'); CdA = CdA.Cdavg; % Coefficients of drag versus mach number

if exist('trade', 'var') % Input handling for trade study
    n = 20;
    lb = [0.004, 1e-5, 0.2, 0.03, 0.02, 3]; % Hardcoded lower bounds
    ub = [0.014, 4e-5, 0.6, 0.08, 0.04, 10]; % Hardcoded upper bounds
    x_pert = ones(n,1).*x;
    x_pert(:,trade) = linspace(lb(trade), ub(trade), n);
    
else % Input handling for design optimization
    x_pert = x./[6,40,33,1,1,3]; % Rescale to unit values MUST BE THE SAME HERE AS IN OPTIMIZATION!!!
    lb = [0.004, 10e-6, 0.2, 0.030, 0.02, 3]; % Lower Bound
    ub = [0.014, 40e-6, 0.6, 0.080, 0.04, 10]; % Upper Bound
    x_pert = (ub-lb).*x_pert+lb; % Scale input to physically meaningful values
end

% Declare additional non-design model input
del_time = 0.01; % (s) Time step
time_max = 30; % (s) Maximum model run time
fill_level = 0.6; % (%) Initial liquid volume fraction in tank
p_tank_init = 5000000; % (Pa) Initial tank pressure
p_feed = 100000; % (Pa) Feed system pressure drop
rho_f = 930; % (kg/m^3) Fuel density
a = 0.000155; % Regression rate constant
n = 0.5; % Regression rate exponent
launchAlt = 1400; % (m) Elevation of Truth or Consequences, New Mexico
T_amb = 301; % (K) 28degC, average of high and low temperatures for Truth or Consequences, New Mexico in June
p_amb = 85600; % (Pa) Pressure at 1400m ASL (Truth or Consequences, New Mexico)
zeta_d = 1.05; % Nozzle discharge correction factor
zeta_cstar = 0.90; % Characteristic velocity correction factor
zeta_CF = 0.90; % Nozzle coefficient correction factor
zeta_Cd = 1; % Drag coefficient variability
p.mass_payload = 4; % (kg) it's a payload
rail_length = 9; % (m) Launch rail length

% Place values into struct for function handling
p.parameters = [del_time, time_max, fill_level, p_tank_init, p_feed, rho_f, a, n, T_amb, p_amb,...
    zeta_d, zeta_cstar, zeta_CF];
p.CEA = CEA;
p.N2Osat = N2Osat;
p.CdA = CdA*zeta_Cd;
p.launchAlt = launchAlt;
p.p_amb = p_amb;


len = size(x_pert);
for k = 1:len(1)
    %% MODEL CALCS %%
    p.design = x_pert(k,:); % Simulate with perturbed design variable

    % Hybrid Motor Model
    try
        p = MuleSim3Opt(p);
    catch ME
        vel_off_rail(k) = NaN; % Error gets thrown if chamber pressure or mass flux is too high, return NaN
        max_alt(k) = NaN;
        I_sp(k) = NaN;
        Thrust_max(k) = NaN;
        pcc_max(k) = NaN;
        continue
    end

    % Structure Estimation
    p.dia = p.d_port_f + 0.05; % Rocket outer diameter is final fuel grain port diameter plus 5cm
    
    if p.dia > 0.1397 % If rocket diameter is larger than skookum-sized
        vel_off_rail(k) = NaN;
        max_alt(k) = NaN;
        I_sp(k) = NaN;
        continue
    end
    
    p.dia = 0.1397; % Set rocket diameter to 5.5"
    p.mass_structural = 0.434866*p.mass_prop + 259.0147*p.dia - 13.2621; % Calculate structural mass from empirically fitted equation
    p.mass = p.mass_structural + p.mass_prop + p.mass_payload;

    % Trajectory Model
    p = suborbitOpt(p);

    % Output
    si = find(p.alt ~= 0, 1, 'first'); % First index that is not zero
    ei = find(p.alt > 9, 1, 'first'); % First index that is greater than nine
    
    vel_off_rail(k) = interp1(p.alt(si:ei), p.vel(si:ei), rail_length); % Calculate off-the-rail velocity
    max_alt(k) = p.apogee; % Maximum altitude (apogee)
    I_sp(k) = 1000 - p.I_sp; % Specific Impulse
    Thrust_max(k) = p.Thrust_max; % Maximum Thrust
    pcc_max(k) = p.pcc_max; % Maximum chamber pressure
    accel_max(k) = max(p.accel); % Maximum mach number

    % If acceleration is too large
    if accel_max(k) > 100 
        I_sp(k) = NaN;
    end
    
    % Penalize objective function if not at target altitude
%     I_sp(k) = I_sp(k)*(1 + (abs(3048 - max_alt(k))/3048)^1.5);
    
    
end

%% PLOTTY PLOT PLOT %%
% 
% figure
% hold on
% plot(x_pert(:,trade), I_sp)
% title('Parametric Study - Specific Impulse')
% xlabel('Nozzle Area Ratio')
% ylabel('Specific Impulse (s)')
% saveas(gcf, ['Output\Parametric_Study\Isp_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
% hold off
% 
% 
% figure
% hold on
% plot(x_pert(:,trade), max_alt)
% title('Parametric Study - Maximum Altitude')
% xlabel('Nozzle Area Ratio')
% ylabel('Altitude (m)')
% saveas(gcf, ['Output\Parametric_Study\alt_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
% hold off
% 
% 
% figure
% hold on
% plot(x_pert(:,trade), pcc_max*1e-6)
% title('Parametric Study - Maximum Pressure')
% xlabel('Nozzle Area Ratio')
% ylabel('Combustion Chamber Pressure (MPa)')
% saveas(gcf, ['Output\Parametric_Study\pcc_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
% hold off
% 
% 
% figure
% hold on
% plot(x_pert(:,trade), Thrust_max)
% title('Parametric Study - Maximum Thrust')
% xlabel('Nozzle Area Ratio')
% ylabel('Thrust (N)')
% saveas(gcf, ['Output\Parametric_Study\Thrust_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
% hold off

t = 0:del_time:(length(p.accel)-1)*del_time;
ylabels = {"Altitude (m)"; "Velocity (m/s)"; "Acceleration (m/s^2)"};
[ax, hlines] = multiplotyyy({t, p.alt},{t, p.vel},{t, p.accel},ylabels);
title('Ramses-1');
xlabel(ax(1),'Time (s)');

end


%                                             &@&.                       
%                                          @@      /@                    
%                               %@@@@@@,  @&    @%   %(                  
%                           (@%         @@@        @                     
%              ,&@@@@@@@@@@@.         @@&         @#                     
%          *@@@@@@&      @/         @@,       ,&,  /@@@.                 
%         @@@@@%        @    &@@@@@@.                 @@%                
%        #@@@@@        @..@*    @@                     @@                
%        *@@@@@        @,    (@/                      &@,                
%         @@@@@@          @@.         *@@@@@,        #@#                 
%          @@@@@@    (@#           #@@      @       @@.                  
%            @@@@@@  .&@@@@@@    @@ @      @/     /@&                    
%             #@@@@@@.    #@   &@  @      @     @@/  #@,                 
%               .@@@@@@@. @@  @@@  @    @.   @@%     @@@%                
%               @  @@@@@@@@@ % @  ,   @%@@@*         #@@@                
%             /#      %@@@@@@@@@.                    @@@@/                       
%            /%           @@@@@@@@@@@@,           (@@@@@@                
%             @          *@.  *@@@@@@@@@@@@@@@@@@@@@@@@@                 
%            @/      .@@            ,&@@@@@@@@@@@@@@@                    
%           @    @@,                                                     
%          @@%                                                           
