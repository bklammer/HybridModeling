function [] = MuleSim2(CEA_f, input_f)
%MULESIM2 UVR Hybrid Rocket Motor Model
% 'input_f' is a string containing the name of a .xlsx file containing miscellaneous motor data
% 'CEA_f' is a string containing the name of a .mat file containing data for combustion reaction and thermodynamic property lookup tables
% All input files must be in the same folder as the 'MuleSim.m' or have the folder path included.

% Benjamin Klammer - 2019
% All calculations performed in metric units

%% INPUTS %%
% Load input files
if nargin == 0 % Set default files to read from in case no input
    input_f = 'MuleSim2INPUT';
    CEA_f = 'MuleSim2CEA';
elseif nargin == 1
    input_f = 'MuleSim2INPUT';
elseif nargin > 2
    error('Invalid number of inputs');
end
input = readtable([input_f, '.xlsx']);
CEA = load([CEA_f, '.mat']);
N2Osat = load('N2Osat.mat');

% Pre-declare variables before they're poofed into existence by eval
[del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    vel_cc_init, d_th, p_atm, A_ratio_nozzle, V_tank_eps, u_tank_eps, ...
    A_ratio_nozzle_eps] = deal(NaN);

for k = 1:length(input.Symbol)
    if ~exist(input.Symbol{k}, 'var')
        error('Variable ''%s'' in the input spreadsheet is not declared in the MATLAB code', input.Symbol{k})
    end
    eval([input.Symbol{k} '= input.Value(' num2str(k) ');']); % Read input data from spreadsheet
end
tic; % Start recording the simulation runtime

% Combustion Parameters
R_const = 8.3144598; % (J/mol*K)
ox_k = fieldnames(CEA.OX); % Read oxidizer info from CEA file
M_ox = zeros(size(ox_k));
for k = 1:length(ox_k)
    M_ox(k) = CEA.OX.(ox_k{k}).M/1000; % (kg/mol)
    Parameter_OX{k} = ['Oxidizer ', num2str(k), ' Mass Fraction'];
    Symbol_OX{k} = ox_k{k};
    Value_OX{k} = CEA.OX.(ox_k{k}).wtfrac;
    Units_OX{k} = '-';
end
fuel_k = fieldnames(CEA.FUEL); % Read fuel info from CEA file
M_f = zeros(size(fuel_k));
h_f = zeros(size(fuel_k));
m_f_wtfrac = zeros(size(fuel_k));
for k = 1:length(fuel_k)
    M_f(k) = CEA.FUEL.(fuel_k{k}).M/1000; % (kg/mol)
    h_f(k) = CEA.FUEL.(fuel_k{k}).h; % (J/mol)
    m_f_wtfrac(k) = CEA.FUEL.(fuel_k{k}).wtfrac; % (kg/kg)
    Parameter_FUEL{k} = ['Fuel ', num2str(k), ' Mass Fraction'];
    Symbol_FUEL{k} = fuel_k{k};
    Value_FUEL{k} = CEA.FUEL.(fuel_k{k}).wtfrac;
    Units_FUEL{k} = '-';
end
products_k = fieldnames(CEA.PRODUCTS); % Read product info from CEA/thermo files
m_k_init = zeros(size(products_k));
M_k = zeros(size(products_k));
for k = 1:length(products_k)
    M_k(k) = CEA.PRODUCTS.(products_k{k}).prop(3)/1000; % (kg/mol)
end
Q_dot_cc = 0; % (W) Ignore the effect of heat losses to the chamber walls (for now)
Q_dot_tank = 0; % (W) Ignore the effect of heat losses to the chamber walls (for now)


%% INITIAL CALCULATIONS %%
% Calculate initial thermodynamic properties from tank pressure
thermo_sat_init = thermoSat(p_tank_init, 'p', {'T', 'rho_liq', 'rho_vap', 'u_liq', 'u_vap'});
thermo_sat_init = num2cell(thermo_sat_init);
[T_tank, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_init{:};
x_tank = (V_tank/m_ox_tank_init - 1/rho_liq)/(1/rho_vap - 1/rho_liq); % Calculate vapour mass fraction
U_tank = m_ox_tank_init*(x_tank*u_vap + (1-x_tank)*u_liq); % (J) Calculate total internal energy

% Calculation of constants
C_inj = n_inj * Cd_inj * (pi()*(d_inj/2)^2); % (m^2) Injector coefficient
A_th = pi()*(d_th/2)^2; % (m^2) Nozzle Throat Area
R_k = (R_const ./ M_k)'; % (J/kg*K)
h_f = sum(h_f.*m_f_wtfrac ./ M_f); % (J/kg)

% Calculate initial mass of oxidizer filling combustion chamber
Vol_init = pi()*L*(d_port_init/2)^2; % (m^3)
R_init = R_const / M_k(strcmp(products_k, 'N2')); % (J/kg*K)
m_init = (p_cc_init*Vol_init) / (R_init*T_cc_init); % (kg) 
m_k_init( strcmp(products_k, 'N2') ) = m_init; % (kg) Set initial contents of combustion chamber to be Nitrogen

% Setting initial conditions
r_cc = d_port_init/2; % (m) Combustion chamber port radius
m_k(1,:) = m_k_init'; % (kg) Array storing mass of each product
T_cc = T_cc_init; % (K) Combustion chamber temperature
vel_cc = vel_cc_init; % (m/s) Combustion chamber velocity
m_ox_tank = m_ox_tank_init; % (kg) Mass of oxidizer in tank
Ma = 3; % Initialize the mach number to three


%% MODEL CALCULATIONS %%
t = 1;
try % Start try block to catch errors that occur in the main loop (for debugging)
while 1 % Run until break
    if any(m_k(t,:) < 0) % Check that there are no negative masses
        error('MuleSim2:negMass', 'One of the product masses is negative')
    end

    % THERMODYNAMIC STATE OF COMBUSTION CHAMBER %
    A(t) = pi()*r_cc(t)^2;
    Vol(t) = A(t)*L;
    m_cc(t) = sum(m_k(t,:));
    rho_cc(t) = m_cc(t)/Vol(t);
    [cp_k(t,:), h_k(t,:)] = thermoProp(products_k, T_cc(t));
    h_cc(t) = sum(h_k(t,:).*m_k(t,:)) / m_cc(t);
    R_cc(t) = sum(R_k.*m_k(t,:)) / m_cc(t);
    cp_cc(t) = sum(cp_k(t,:).*m_k(t,:)) / m_cc(t);
    cv_cc(t) = cp_cc(t) - R_cc(t);
    k_cc(t) = cp_cc(t)/cv_cc(t);
    u_cc(t) = cv_cc(t)*T_cc(t);
    p_cc(t) = rho_cc(t)*R_cc(t)*T_cc(t);
    
    
    if x_tank(t) < 1
        % OXIDIZER TANK CALCULATIONS FOR SATURATED LIQUID AND VAPOUR %
        if abs(Verror(T_tank(t))) > V_tank_eps
            T_tank(t) = secant(@Verror, T_tank(t));
        end
        
        % Calculation of tank thermo properties
        thermo_sat_tank = thermoSat(T_tank(t), 'T', {'p', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap', 'u_liq', 'u_vap'});
        thermo_sat_tank = num2cell(thermo_sat_tank);  
        [p_tank(t), h_liq, h_vap, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_tank{:};
        x_tank(t) = (U_tank(t)/m_ox_tank(t) - u_liq)/(u_vap - u_liq);
        u_tank(t) = x_tank(t)*u_vap + (1-x_tank(t))*u_liq;
        h_tank(t) = x_tank(t)*h_vap + (1-x_tank(t))*h_liq;
        rho_tank(t) = 1/(x_tank(t)/rho_vap + (1-x_tank(t))/rho_liq);
        h_ox_in(t) = h_liq + 1.38548e+06; % Convert from NIST enthalpy convention to CEA enthalpy convention 
        
        % Thermo properties downstream of injector for HEM model
        s_liq = thermoSat(T_tank(t), 'T', {'s_liq'});
        thermo_sat_inj = thermoSat(p_cc(t), 'p', {'s_liq', 's_vap', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap'});
        thermo_sat_inj = num2cell(thermo_sat_inj);
        [s_liq_2, s_vap_2, h_liq_2, h_vap_2, rho_liq_2, rho_vap_2] = thermo_sat_inj{:};
        x_2(t) = (s_liq - s_liq_2)/(s_vap_2 - s_liq_2); % Mass fraction downstream assuming isentropic expansion at saturation (only use for HEM mass flow)
        rho_ox_2(t) = 1/(x_2(t)/rho_vap_2 + (1 - x_2(t))/rho_liq_2);
        h_2(t) = x_2(t)*h_vap_2 + (1 - x_2(t))*h_liq_2;

        % Calculate mass flow through the injector
        G_SPI(t) = sqrt(2*rho_liq*(p_tank(t)-p_cc(t)));
        G_HEM(t) = rho_ox_2(t)*sqrt(2*(h_liq-h_2(t)));
        m_dot_ox_in(t) = C_inj * (G_SPI(t) + G_HEM(t))/2;
        
        del_m_ox_tank = -m_dot_ox_in(t)*del_time;
        del_U_tank = (-m_dot_ox_in(t)*h_liq + Q_dot_tank)*del_time;
        
        
    else
        % OXIDIZER TANK CALCULATIONS FOR VAPOUR ONLY %
        rho_tank(t) = m_ox_tank(t)./V_tank;
        u_tank(t) = U_tank(t)/m_ox_tank(t);
        if abs(uerror(T_tank(t))) > u_tank_eps
            T_tank(t) = secant(@uerror, T_tank(t));
        end
        
        thermo_span_tank = thermoSpanWagner({'p','h'}, rho_tank(t), T_tank(t));
        p_tank(t) = thermo_span_tank(1);
        h_tank(t) = thermo_span_tank(2) + 7.3397e+05; % Convert from Span-Wagner enthalpy convention to NIST
        h_ox_in(t) = h_tank(t) + 1.38548e+06; % Convert from NIST enthalpy convention to CEA
        m_dot_ox_in(t) = C_inj*sqrt(2*rho_tank(t)*(p_tank(t)-p_cc(t))); % Incompressible fluid assumption (better than nothing)
        
        del_m_ox_tank = -m_dot_ox_in(t)*del_time;
        del_U_tank = (-m_dot_ox_in(t)*h_tank(t) + Q_dot_tank)*del_time;
    end
    
    if p_cc(t) > 0.9*p_tank(t) % Check that the chamber pressure is sufficiently lower than the tank pressure
        error('MuleSim2:highPressure', 'The combustion chamber pressure is too high (p_cc=%7.0f Pa) \n', p_cc(t))
    end
    
    % NOZZLE AND THRUST CALCULATIONS %
    T_stag(t) = T_cc(t) + vel_cc(t).^2/(2.*cp_cc(t)); % (K) Stagnation temperature is calculated since flow has non-negligible velocity
    p_stag(t) = p_cc(t)*(T_stag(t)/T_cc(t)).^(k_cc(t)/(k_cc(t)-1)); % (Pa) Stagnation pressure is calculated since flow has non-negligible velocity
    p_crit = p_atm/(2/(k_cc(t)+1))^(k_cc(t)/(k_cc(t)-1)); % Critical pressure required for choked flow
    if p_stag(t) > p_crit % As long as the choked flow assumption is valid
        m_dot_out(t) = A_th.*p_stag(t)./sqrt(T_stag(t))*sqrt(k_cc(t)./R_cc(t)).*(2./(k_cc(t)+1)).^((k_cc(t)+1)./(2*k_cc(t)-2));
    else
        rho_stag(t) = p_stag(t)/(R_cc(t)*T_stag(t));
        m_dot_out(t) = A_th*sqrt(2*rho_stag(t)*(p_stag(t)-p_atm)); % Incompressible fluid mass flow
    end
    
    vel_out(t) = m_dot_out(t) / (rho_cc(t)*A(t)); % Velocity from exit mass flow
    if vel_out(t) > L/del_time % Check that mass can't pass through the chamber without the model seeing it
        error('MuleSim2:highVel', 'The outlet velocity (vel_out = %3.0f m/s) is too high, try decreasing ''del_time''', vel_out(t))
    end
    
    if abs(Aerror(Ma(t))) > A_ratio_nozzle_eps % Nozzle mach number solver
        Ma(t) = secant(@Aerror, Ma(t));
    end
    p_exit(t) = p_stag(t)./(1+(k_cc(t)-1)./2.*Ma(t).^2).^(k_cc(t)/(k_cc(t)-1)); % (Pa) Pressure of gas at exit of nozzle
    if p_exit(t) < 0.3*p_atm %x_tank(t) >= 1 % Physically very wrong, but the condition used to hide the transition that results from overexpansion
        p_exit(t) = p_stag(t)*(2/(k_cc(t)+1))^(k_cc(t)/(k_cc(t)-1)); % Throat pressure
        vel_exit(t) = sqrt(2*k_cc(t)*R_cc(t)*T_stag(t)/(k_cc(t)+1)); % Throat velocity
        A_ratio_nozzle_eff = 1;
    else
        T_exit(t) = T_stag(t)./(1+(k_cc(t)-1)./2.*Ma(t).^2); % (K) Temperature of gas at exit of nozzle
        vel_exit(t) = Ma(t).*sqrt(k_cc(t).*R_cc(t).*T_exit(t)); % (m/s) Velocity of gas at exit of nozzle
        A_ratio_nozzle_eff = A_ratio_nozzle; % Effective nozzle area ratio is actual nozzle area ratio
     end
    F_thrust(t) = m_dot_out(t).*vel_exit(t) + (p_exit(t)-p_atm).*A_th.*A_ratio_nozzle_eff; % (N) Rocket Motor Thrust!!!
    
    
    % COMBUSTION CHAMBER CALCULATIONS %
    vel_cc(t) = vel_out(t); % Set cc velocity to nozzle inlet velocity
    if t > 1 % As long as not the first time step
        del_vel_cc(t) = vel_cc(t)-vel_cc(t-1); % Use backwards difference to calculate change in velocity
    else
        del_vel_cc(t) = 0; % Set to zero if first time step
    end
    
    % Regression
    del_r_cc(t) = a*(m_dot_ox_in(t) / A(t))^n * del_time; % Regression rate law
    m_dot_f(t) = 2*pi()*r_cc(t)*L * rho_f * (del_r_cc(t)/del_time); % Fuel mass flow rate from geometry
    
    % Mass balance
    m_dot_i(t) = m_dot_f(t) + m_dot_ox_in(t); % Total mass flow in
    OF(t) = m_dot_ox_in(t)./m_dot_f(t); % Calculate Oxidizer-Fuel Ratio
    alpha_k(t,:) = massFrac(OF(t), p_cc(t)); % Interpolate product mass fractions from CEA lookup table
    del_m_k(t,:) = (alpha_k(t,:).*m_dot_i(t) - m_k(t,:).*vel_cc(t)/L) .* del_time; % Product mass balances (inflow from combustion, outflow to nozzle)
    del_m_cc(t) = sum(del_m_k(t,:)); % Total mass balance
     
    % Energy balance
    del_U_cc(t) = (m_dot_ox_in(t)*h_ox_in(t) + m_dot_f(t)*h_f - m_dot_out(t)*(h_cc(t) + (vel_out(t)^2)/2) + Q_dot_cc) * del_time; % The energy equation (adds constant for h_ox_in to account for different standards of enthalpy)
    del_T_cc(t) = (del_U_cc(t) - (u_cc(t) + (vel_cc(t)^2)/2)*del_m_cc(t) - m_cc(t)*vel_cc(t)*del_vel_cc(t)) / (m_cc(t)*cv_cc(t));
    
    
    % ITERATE FORWARD IN TIME %
    % If less than max burn time and more than 3% of oxidizer is left, step forward in time
    if t < time_max/del_time && m_ox_tank(t)/m_ox_tank_init > 0.03
        % Oxidizer Tank Values
        m_ox_tank(t+1) = m_ox_tank(t) + del_m_ox_tank;
        U_tank(t+1) = U_tank(t) + del_U_tank;
        x_tank(t+1) = x_tank(t);
        T_tank(t+1) = T_tank(t);
        % Combustion Chamber Values
        r_cc(t+1) = r_cc(t) + del_r_cc(t);
        m_k(t+1,:) = m_k(t,:) + del_m_k(t,:);
        vel_cc(t+1) = vel_cc(t) + del_vel_cc(t);
        T_cc(t+1) = T_cc(t) + del_T_cc(t);
        Ma(t+1) = Ma(t);
        t = t+1;
    else
        break % Exit loop if no more time or oxidizer
    end
end
catch ME
    fprintf('Something went wrong at time t=%f \n', t*del_time);
    rethrow(ME)
end
time = 0:del_time:(del_time*(t-1)); % Create a time vector


%% OUTPUTS %%
% Major output calculations
r_dot = del_r_cc./del_time; % (m/s)
d_port_f = 2*r_cc(end); % (m)
m_out = trapz(m_dot_out)*del_time; % (kg)
m_f = trapz(m_dot_f)*del_time; % (kg)
m_ox = m_out-m_f; % (kg)
I_tot = trapz(F_thrust)*del_time; % (N*s)
v_e = I_tot/m_out; % (m/s)
I_sp = v_e/9.81; % (s)
c_star = mean(p_cc)*A_th/mean(m_dot_out); % (m/s)

% Save metadata of simulation
time_sim = toc; % Save how long the simulation took (minus loading data and plotting figures)
time_stamp = datestr(datetime); % Save the current time and date
computer_name = computer; % Save which computer the code is running on
function_name = mfilename; % Save the name of the function that is being called


%% PLOTS %%
figure(1)
hold on
plot(time, F_thrust./4.44822)
plot(time, p_tank./6894.76)
plot(time, p_cc./6894.76)
title('Boundless - University of Washington Results')
xlabel('Time (s)')
ylabel('Thrust (lbf) and Pressure (psi)')
legend('Thrust', 'Tank Pressure', 'Chamber Pressure')
axis([-1 20 -50 1400])
hold off

% figure(2)
% hold on
% plot(time, r_cc*1000)
% title('Fuel Grain Radius vs. Time')
% xlabel('Time (s)')
% ylabel('Radius (mm)')
% axis([0 inf 0 inf])
% hold off
% 
% figure(3)
% hold on
% plot(time, T_cc)
% title('Combustion Chamber Temperature vs. Time')
% xlabel('Time (s)')
% ylabel('Temperature (K)')
% hold off
% 
% figure(4)
% hold on
% plot(time, F_thrust)
% title('Rocket Thrust vs. Time')
% xlabel('Time (s)')
% ylabel('Force (N)')
% hold off
% 
% figure(3)
% hold on
% plot(time, m_dot_out)
% title('Total Mass Flow vs. Time')
% xlabel('Time (s)')
% ylabel('Mass Flow (kg/s)')
% hold off
% 
% figure(4)
% hold on
% plot(time, r_dot*1000)
% title('Fuel Grain Regression Rate vs. Time')
% xlabel('Time (s)')
% ylabel('Regression Rate (mm/s)')
% axis([0 inf 0 inf])
% hold off

%% SAVE OUTPUT TO SPREADSHEETS %%
% Add output data to table
Parameter_out = {'Burn Time'; 'Peak Thrust'; 'Average Thrust'; ...
    'Total Impulse'; 'Specific Impulse'; 'Average O/F Ratio'; ...
    'Average Oxidizer Mass Flow'; 'Average Fuel Mass Flow'; ...
    'Total Oxidizer Mass Consumed'; 'Total Fuel Mass Consumed'; ...
    'Average Regression Rate'; 'Final Port Diameter'; ...
    'Maximum Combustion Chamber Pressure'; ...
    'Maximum Combustion Chamber Temperature'}; 
Symbol_out = {'time(end)'; 'max(F_thrust)'; 'mean(F_thrust)'; 'I_tot'; ...
    'I_sp'; 'mean(OF)'; 'mean(m_dot_ox_in)'; 'mean(m_dot_f)'; 'm_ox'; ...
    'm_f'; 'mean(r_dot)'; 'd_port_f'; 'max(T_cc)'; 'max(p_cc)'};
Value_out = {time(end); max(F_thrust); mean(F_thrust); I_tot; I_sp; ...
    mean(OF); mean(m_dot_ox_in); mean(m_dot_f); m_ox; m_f; mean(r_dot); ...
    d_port_f; max(T_cc); max(p_cc)};
Units_out = {'s'; 'N'; 'N'; 'N*s'; 's'; '-'; 'kg/s'; 'kg/s'; 'kg'; ...
    'kg'; 'm/s'; 'm'; 'K'; 'Pa'};

% Add metadata to table
Parameter_meta = {'MATLAB Function Name'; 'Date and Time Simulated'; ...
    'CEA File Name'; 'Input File Name'; 'Elapsed Simulation Time'; 'Operating System and Version'};
Symbol_meta = {'function_name'; 'time_stamp'; 'CEA_f'; 'input_f'; 'time_sim'; 'computer_name'};
Value_meta = {function_name; time_stamp; CEA_f; input_f; time_sim; computer_name};
Units_meta = {'-'; '-'; '-'; '-'; 's'; '-'};

% Output to table
Parameter = [Parameter_meta; Parameter_OX'; Parameter_FUEL'; input.Parameter; Parameter_out];
Symbol = [Symbol_meta; Symbol_OX'; Symbol_FUEL'; input.Symbol; Symbol_out];
Value = [Value_meta; Value_OX'; Value_FUEL'; num2cell(input.Value); Value_out];
Units = [Units_meta; Units_OX'; Units_FUEL'; input.Units; Units_out];

% Combine data to table and write to file
output = table(Parameter, Symbol, Value, Units);
filename = ['Output\', CEA_f, '_summary_', datestr(now, 'yyyy-mm-dd'), '.csv']; % Create a filename for the table based on the function name and the date, save it to the 'Output' folder
writetable(output, filename);

% Output Time-series data
output_time = table(['[s]'; num2cell(time')], ['[N]'; num2cell(F_thrust')], ...
    ['[K]'; num2cell(T_cc')], ['[Pa]'; num2cell(p_cc')], ['[-]'; num2cell(OF')], ...
    ['[K]'; num2cell(T_tank')], ['[Pa]'; num2cell(p_tank')], ...
    ['[-]'; num2cell(x_tank')]);
output_time.Properties.VariableNames = {'Time', 'Thrust', ...
    'T_cc', 'p_cc', 'OF', 'T_tank', 'p_tank', 'x_tank'};
filename = ['Output\', CEA_f, '_time_', datestr(now, 'yyyy-mm-dd'), '.csv']; % Create a filename for the table based on the function name and the date, save it to the 'Output' folder
writetable(output_time, filename);

%% SUBFUNCTIONS %%
    function [alpha] = massFrac(OF, p)
    %MASSFRAC Calculates the mass fractions of products for a given reaction
    %   'OF' is the oxidizer to fuel weight ratio
    %   'p' is the pressure of the combustion chamber
    %   'alpha' is an array of weight fractions
        p = p/100000; % Convert Pa to bar for lookup table
        if OF < CEA.OF(1) % Check that query is inside bounds
            error('massFrac:lowOF', 'Outside O/F range: \n OF = %f', OF)
        elseif OF > CEA.OF(end)
            error('massFrac:highOF', 'Outside O/F range: \n OF = %f', OF)
        elseif p < CEA.P(1)
            error('massFrac:lowP', 'Outside pressure range: \n p = %f bar', p)
        elseif p > CEA.P(end)
            error('massFrac:highP', 'Outside pressure range: \n p = %f bar', p)
        end
        for mm = 1:length(CEA.OF)
            if OF < CEA.OF(mm) % Find first index that is larger than OF
                break
            end
        end
        for nn = 1:length(CEA.P)
            if p < CEA.P(nn) % Find first index that is larger than p
                break
            end
        end
        names = fieldnames(CEA.PRODUCTS);
        alpha = zeros(size(names));
        for kk = 1:length(names) % Iterate through all products
            % Do some bilinear interpolation (equations from wikipedia)
            OF_int = [CEA.OF(mm)-OF, OF-CEA.OF(mm-1)];
            P_int = [CEA.P(nn)-p; p-CEA.P(nn-1)];
            C_int = 1/((CEA.OF(mm)-CEA.OF(mm-1))*(CEA.P(nn)-CEA.P(nn-1)));
            alpha_int = CEA.PRODUCTS.(names{kk}).MASSFRAC(mm-1:mm,nn-1:nn);
            alpha(kk) = C_int.* OF_int*alpha_int*P_int; % Mass fraction of product k
        end
        alpha = alpha'; % Transpose alpha so that it fits the same format in the main function
    end


    function [cp, h] = thermoProp(products, T)
    %THERMOPROP Calculates the thermodynamic properties of given set of species
    %   'products' is a cell array of the names of the product species
    %   'T' is the temperature at which the properties are to be evaluated (K)
    %   'cp' is the specific heat capacity at constant pressure (J/kg*K)
    %   'h' is the specific enthalpy (J/kg)
        cp = zeros(size(products)); % Initialize empty arrays
        h = zeros(size(products));
        for kk = 1:length(products)
            T_temp = T;
            tp = CEA.PRODUCTS.(products{kk}); % Save product field as local variable for code readibility
            M = tp.prop(3)/1000; % (kg/mol)
            Rk = R_const/M; % (J/kg*K)
            if T_temp < tp.T(1)-1 % Error if outside of temperature range
                error('thermoProp:lowTemp', 'Out of temperature range: \n T = %5.0f K ', T)
                T_temp = tp.T(1);
            elseif T_temp > tp.T(end)+1
                error('thermoProp:highTemp', 'Out of temperature range: \n T = %5.0f K ', T)
                T_temp = tp.T(end)-1;
            end
            for mm = 1:tp.prop(1) % Find appropriate temperature range
                if T_temp < tp.T(mm,2)
                    break
                end
            end
            cp(kk) = Rk*sum( tp.a(mm,:) .* (T_temp .^ (tp.n(mm,:))) ); % (J/kg*K) Array-wise calculation of cp, see NASA's CEA manual for details
            h_temp = (tp.a(mm,:)./(tp.n(mm,:)+1)) .* (T_temp .^ (tp.n(mm,:)+1)); % Uses the formula for integration of a polynomial, since delh = integral(cp)dT
            h_temp(~isfinite(h_temp)) = tp.a(mm, ~isfinite(h_temp)).* log(T_temp); % Replaces the infinite term with ln(T) since polynomial integration formula doesn't work for 1/T
            h(kk) = Rk*(sum(h_temp) + tp.b(mm,1)); % (J/kg)
        end
    end


    function [varO] = thermoSat(varI, varIname, varOname)
    %THERMOSAT returns thermodynamic properties at saturation
    %   'varI' is a double specifying the value of thermodynamic property inputted
    %   'varIname' is a string specifying the given input thermodynamic property 
    %   'varOname' as a cell array specifying the desired output thermodynamic property
        varO = zeros(size(varOname));
        for mm = 1:length(varOname) % For each desired output
            % Locate which variable column is input, and which is output
            ii = find( (varIname == N2Osat.meta(1,:)), 1, 'first');
            jj = find( (varOname(mm) == N2Osat.meta(1,:)), 1, 'first');
            % Interpolate between the two values
            for kk = 1:length(N2Osat.data(:,ii))
                if varI < N2Osat.data(kk,ii)
                    break
                end
            end
            varO(mm) = (varI - N2Osat.data(kk-1,ii))./(N2Osat.data(kk,ii)-N2Osat.data(kk-1,ii)).*(N2Osat.data(kk,jj)-N2Osat.data(kk-1,jj)) + N2Osat.data(kk-1,jj);
        end
    end


    function [x] = secant(fun,x1)
    %SECANT is a zero-finding function based on the secant method
    %   'fun' is the function handle for which the zero is desired
    %   'x1' is the initial guess
        x_eps = x1*0.005; % Set the tolerance to be 0.5% of initial guess
        x2 = x1-x1*0.01; % Set a second point 1% away from the original guess
        F1 = fun(x1); % Evaluate function at x1
        F2 = fun(x2); % Evaluate function at x2
        kk = 1; % Set up counter
        kk_max = 1000;
        while abs(x2-x1)>=x_eps && kk<kk_max % While error is too large and counter is less than max
            x3 = x2 - F2*(x2-x1)/(F2-F1);
            x1 = x2; % Move everything forward
            x2 = x3;
            F1 = F2;
            F2 = fun(x2);
            kk = kk+1;
        end
        x = x2;
    end


    % SET UP FUNCTIONS TO BE ITERATIVELY SOLVED %
    function V = Verror(T) % Finds the difference between the estimated and actual tank volume
        thermo = thermoSat(T, 'T', {'rho_liq','rho_vap','u_liq','u_vap'});
        thermo = num2cell(thermo);
        [rho_l,rho_v,u_l,u_v] = thermo{:};
        x = (U_tank(t)/m_ox_tank(t) - u_l)/(u_v - u_l);
        V = m_ox_tank(t)*((1-x)/rho_l + x/rho_v) - V_tank;
    end
    
    function U = uerror(T) % Finds the difference between the estimated and actual tank internal energy
        U = (thermoSpanWagner({'u'}, rho_tank(t), T) + 7.3397e+5) - u_tank(t);
    end

    function A = Aerror(M) % Finds the difference between the estimated and actual nozzle ratio
    	A = (1/M^2)*(2./(k_cc(t)+1).*(1+(k_cc(t)-1)./2.*M.^2)).^((k_cc(t)+1)./(k_cc(t)-1)) - A_ratio_nozzle^2;
    end

end
