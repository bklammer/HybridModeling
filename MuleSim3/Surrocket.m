function [obj, state] = Surrocket(x)
%SURROCKET Models a hybrid rocket from design to apogee
% For Design Optimization
%   'x' is a design vector containing the following variables:
%       'V_tank' - oxidizer tank volume (m^3)
%       'C_inj' - effective injector area (m^2)
%       'L' - fuel grain length (m)
%       'd_port_init' - initial fuel grain port diameter (m)
%       'd_th' - nozzle throat diameter (m)
%       'A_ratio_nozzle' - nozzle area ratio as a fraction
%   'notscaled' is a dummy variable to indicate that the input parameters
%   are not scaled
%
% For Uncertainty Analysis
%   'x' is an uncertainty vector containing the following variables:
%       'V_tank' - oxidizer tank volume (m^3)
%       'fill_level' - volume fraction of liquid oxidizer in tank
%       'p_tank_init' - initial oxidizer tank pressure (Pa)
%       'C_inj' - effective injector area (m^2)
%       'L' - fuel grain length (m)
%       'd_port_init' - initial fuel grain port diameter (m)
%       'd_th' - nozzle throat diameter (m)
%       'A_ratio_nozzle' - nozzle area ratio as a fraction
%       'p_feed' - feed system pressure drop (Pa)
%       'rho_f' - solid fuel density (kg/m^3)
%       'a' - regression rate coefficient (m/s)
%       'n' - regression rate exponent
%       'launchAlt' - altitude of launch site (m)
%       'p_amb' - ambient pressure at launch site (Pa)
%       'T_amb' - ambient temperature at launch site (K)
%       'zeta_d' - discharge correction factor
%       'zeta_cstar' - characteristic velocity correction factor
%       'zeta_CF' - thrust coefficient correction factor
%       'C_D' - drag coefficient variability

%% INPUTS %%
% Load input files
state.CEA = load('MuleSim3CEA.mat'); % Combustion product lookup table
state.N2Osat = load('N2Osat.mat'); % Nitrous oxide saturation properties
CdA = load('Cdavg.mat'); CdA = CdA.Cdavg; % Coefficients of drag versus mach number

if length(x)==6 % If used for design optimization
    % Input un-scaling
    x = x./[6,40,33,1,1,3]; % Rescale to unit values MUST BE THE SAME HERE AS IN "OPTIMIZATION.M"!!!
    lb = [0.004, 10e-6, 0.2, 0.030, 0.02, 3]; % Lower Bound
    ub = [0.014, 40e-6, 0.6, 0.080, 0.04, 10]; % Upper Bound
    x = (ub-lb).*x+lb; % Scale input to physically meaningful values
    
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

    % Place values into struct for function handling
    state.design = x;
    state.parameters = [del_time, time_max, fill_level, p_tank_init, ...
        p_feed, rho_f, a, n, T_amb, p_amb, zeta_d, zeta_cstar, zeta_CF];
    state.CdA = CdA*zeta_Cd;
    state.launchAlt = launchAlt;
    state.p_amb = p_amb;
    
elseif length(x) == 7 % If used for Parametric Study
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

    % Place values into struct for function handling
    state.parametricstudy = x(1:6);
    state.parameters = [del_time, time_max, fill_level, p_tank_init, ...
        p_feed, rho_f, a, n, T_amb, p_amb, zeta_d, zeta_cstar, zeta_CF];
    state.CdA = CdA*zeta_Cd;
    state.launchAlt = launchAlt;
    state.p_amb = p_amb;
    
elseif length(x) == 19 % If used for Sensitivity Analysis
    % Place values into struct for function handling
    state.uncertainty = x;
    state.CdA = CdA*x(19); % Adjust drag coefficient
    state.launchAlt = x(13); % Set launch altitude
    state.p_amb = x(14); % Set ambient pressure
end

state.mass_payload = 4; % (kg) it's a payload
state.rail_length = 9; % (m) Launch rail length

%% MODEL CALCULATIONS %%

% Hybrid Motor Model
try
    state = MuleSim3(state);
catch ME
    state.ME = ME; % If error gets thrown, return NaN
    state.I_sp = NaN;
    state.apogee = NaN;
    state.pcc_max = NaN;
    state.Thrust_max = NaN;
    obj = NaN;
    return
end

% Structure Estimation
state.dia = state.d_port_f + 0.05; % Rocket outer diameter is final fuel grain port diameter plus 5cm

if isfield(state, 'design') % Enforce constraint if design is being optimized
    if state.dia > 0.1397 % If rocket diameter is larger than skookum-sized
        state.I_sp = NaN;
        state.apogee = NaN;
        state.pcc_max = NaN;
        state.Thrust_max = NaN;
        obj = NaN;
        return
    end
end

state.dia = 0.1397; % Set rocket diameter to 5.5"
state.mass_structural = 0.434866*state.mass_prop + 259.0147*state.dia - 13.2621; % Calculate structural mass from empirically fitted equation
state.mass = state.mass_structural + state.mass_prop + state.mass_payload;

% Trajectory Model
state = suborbitOpt(state);
si = find(state.alt ~= 0, 1, 'first'); % First index that is not zero
ei = find(state.alt > 9, 1, 'first'); % First index that is greater than nine
state.vel_off_rail = interp1(state.alt(si:ei), state.vel(si:ei), state.rail_length); % Calculate off-the-rail velocity

if isfield(state, 'design') % Enforce constraint if design is being optimized
    % Calculate Objective Function
    obj = 1000 - state.I_sp;
    
    % If acceleration is too large, constraint is violated
    if max(state.accel) > 100 
        obj = NaN;
    end
    % Penalize objective function if not at target altitude
    obj = obj*(1 + (abs(3048 - state.apogee)/3048)^1.5);
    
else
    % If used for uncertainty analysis, return variable of interest
    obj = state.apogee;
end


%% PLOTTY PLOT PLOT %%

% t = 0:del_time:(length(p.accel)-1)*del_time;
% ylabels = {"Altitude (m)"; "Velocity (m/s)"; "Acceleration (m/s^2)"};
% [ax, hlines] = multiplotyyy({t, p.alt},{t, p.vel},{t, p.accel},ylabels);
% title('Ramses-1');
% xlabel(ax(1),'Time (s)');

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
