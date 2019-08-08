function [I_sp, max_alt, vel_off_rail, Thrust_max, pcc_max, p] = SurrocketSensitivity(x)
%SURROCKET Models a hybrid rocket from design to apogee
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
CEA = load('MuleSim3CEA.mat'); % Combustion product lookup table
N2Osat = load('N2Osat.mat'); % Nitrous oxide saturation properties
CdA = load('Cdavg.mat'); CdA = CdA.Cdavg; % Coefficients of drag versus mach number

% Place values into struct for function handling
p.CEA = CEA;
p.N2Osat = N2Osat;
p.uncertainty = x;
p.CdA = CdA*x(19); % Adjust drag coefficient
p.launchAlt = x(13); % Set launch altitude
p.p_amb = x(14);

%% MODEL CALCS %%
% Hybrid Motor Model
try
    p = MuleSim3(p);
catch ME
    vel_off_rail = NaN; % Error gets thrown if combustion pressure is too high, return NaN
    max_alt = NaN;
    I_sp = NaN;
    Thrust_max = NaN;
    pcc_max = NaN;
    return
end

% Structure Estimation
p.dia = p.d_port_f + 0.05; % Rocket outer diameter is final fuel grain port diameter plus 5cm
p.dia = 0.1397; % Set rocket diameter to 5.5"
p.mass_structural = 0.434866*p.mass_prop + 259.0147*p.dia - 13.2621;
p.mass = p.mass_structural + p.mass_prop;

% Trajectory Model
p = suborbitOpt(p);

% Output
si = find(p.alt ~= 0, 1, 'first'); % First index that is not zero
ei = find(p.alt > 9, 1, 'first'); % First index that is greater than nine

vel_off_rail = interp1(p.alt(si:ei), p.vel(si:ei), 9); % Calculate off-the-rail velocity
max_alt = p.apogee; % Maximum altitude (apogee)
I_sp = p.I_sp; % Specific Impulse
Thrust_max = p.Thrust_max;
pcc_max = p.pcc_max;

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
