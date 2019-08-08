function [I_tot, T_curve] = MuleSim3(CEA_f, input_f)
%MULESIM3 UVR Hybrid Rocket Motor Model
% 'input_f' is a string containing the name of a .xlsx file containing miscellaneous motor data
% 'CEA_f' is a string containing the name of a .mat file containing data for combustion reaction
% All input files must be in the same folder as the 'MuleSim3.m' or have the folder path included

% Benjamin Klammer - 2019
% All calculations performed in metric units

%% INPUTS %%
% Load input files
if nargin == 0 % Set default files to read from in case no input
    input_f = 'MuleSim3INPUT';
    CEA_f = 'Phoenix1A_CEA'; % CHANGE ME
elseif nargin == 1
    input_f = 'MuleSim3INPUT';
elseif nargin > 2
    error('Invalid number of inputs');
end
input = readtable([input_f, '.xlsx']);
CEA = load([CEA_f, '.mat']);
N2Osat = load('N2Osat.mat');
valid = load('Motor_Validation.mat');

% Pre-declare variables before they're poofed into existence by eval
[del_time, time_max, V_tank, p_tank_init, p_feed, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, d_f, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    zeta_d, zeta_cstar, zeta_CF, d_th, p_atm, A_ratio_nozzle, graph, save] = deal(NaN);

for k = 1:length(input.Symbol)
    if ~exist(input.Symbol{k}, 'var')
        error('Variable ''%s'' in the input spreadsheet is not declared in the MATLAB code', input.Symbol{k})
    end
    eval([input.Symbol{k} '= input.Value(' num2str(k) ');']); % Read input data from spreadsheet
end

tic; % Start recording the simulation runtime

%% INITIAL CALCULATIONS %%
% Calculate initial thermodynamic properties from tank pressure
thermo_sat_init = thermoSat(p_tank_init, 'p', {'T', 'rho_liq', 'rho_vap', 'u_liq', 'u_vap'});
thermo_sat_init = num2cell(thermo_sat_init);
[T_tank, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_init{:};
x_tank = (V_tank/m_ox_tank_init - 1/rho_liq)/(1/rho_vap - 1/rho_liq); % Calculate vapour mass fraction
u_tank = (x_tank*u_vap + (1-x_tank)*u_liq); % (J/kg*K) Calculate specific internal energy
U_tank = m_ox_tank_init*u_tank; % (J) Calculate total internal energy

% Calculation of constants
R_const = 8.3144598; % (J/mol*K) Gas constant
Q_dot_tank = 0; % (W) Ignore the effect of heat losses to the tank walls (for now)
C_inj = n_inj * Cd_inj * (pi()*(d_inj/2)^2); % (m^2) Injector coefficient
A_th = pi()*(d_th/2)^2; % (m^2) Nozzle Throat Area
V_tank_eps = 0.001*V_tank; % (m^3) Set acceptable tank volume error to 0.1% of tank volume
u_tank_eps = 0.001*u_tank; % (m^3) Set acceptable tank specific internal energy error to 0.1% of initial tank specific internal energy
A_ratio_nozzle_eps = 0.001*A_ratio_nozzle; % (m^3) Set acceptable nozzle area ratio error to 0.1% of nozzle area ratio
m_f_init = rho_f*L*pi()/4*(d_f^2 - d_port_init^2); % (kg) Calculate initial fuel mass from grain geometry

% Setting initial conditions
m_tot = m_f_init + m_ox_tank_init;
m_dot_ox_in = 1; % Initialize 
m_ox_tank = m_ox_tank_init; % (kg) Mass of oxidizer in tank
r_cc = d_port_init/2; % (m) Combustion chamber port radius
p_cc = p_cc_init; % (kg) Array storing mass of each product
p_stag = p_cc_init;
T_cc = T_cc_init; % (K) Combustion chamber temperature
T_stag = T_cc;
R_cc = R_const/0.02897; % (J/kg*K) Set initial fluid in chamber to be air
k_cc = 1.4; % Values for air at 298K
Ma = 3; % Initialize the mach number to three


%% MODEL CALCULATIONS %%
t = 1;
try % Start try block to catch errors that occur in the main loop (for debugging)
while 1 % Run until break
    
    if x_tank(t) < 1
        % OXIDIZER TANK CALCULATIONS FOR SATURATED LIQUID AND VAPOUR %
        Vin.U = U_tank(t);
        Vin.m = m_ox_tank(t);
        Vin.V = V_tank;
        if abs(Verror(T_tank(t), Vin)) > V_tank_eps
            T_tank(t) = secant(@Verror, T_tank(t), Vin);
        end
        
        % Calculation of tank thermo properties
        thermo_sat_tank = thermoSat(T_tank(t), 'T', {'p', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap', 'u_liq', 'u_vap'});
        thermo_sat_tank = num2cell(thermo_sat_tank);  
        [p_tank(t), h_liq, h_vap, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_tank{:};
        x_tank(t) = (U_tank(t)/m_ox_tank(t) - u_liq)/(u_vap - u_liq);
        u_tank(t) = x_tank(t)*u_vap + (1-x_tank(t))*u_liq;
        h_tank(t) = x_tank(t)*h_vap + (1-x_tank(t))*h_liq;
        rho_tank(t) = 1/(x_tank(t)/rho_vap + (1-x_tank(t))/rho_liq);
        rho_discharge = rho_liq;
        h_discharge = h_liq;
        
        % Thermo properties downstream of injector for HEM model
        s_liq = thermoSat(T_tank(t), 'T', {'s_liq'});
        thermo_sat_inj = thermoSat(p_cc(t), 'p', {'s_liq', 's_vap', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap'});
        thermo_sat_inj = num2cell(thermo_sat_inj);
        [s_liq_2, s_vap_2, h_liq_2, h_vap_2, rho_liq_2, rho_vap_2] = thermo_sat_inj{:};
        x_2(t) = (s_liq - s_liq_2)/(s_vap_2 - s_liq_2); % Mass fraction downstream assuming isentropic expansion at saturation (only use for HEM mass flow)
        rho_ox_2(t) = 1/(x_2(t)/rho_vap_2 + (1 - x_2(t))/rho_liq_2);
        h_2(t) = x_2(t)*h_vap_2 + (1 - x_2(t))*h_liq_2;

        % Calculate mass flow through the injector
        G_SPI(t) = sqrt(2*rho_liq*(p_tank(t)-p_cc(t)-p_feed));
        G_HEM(t) = rho_ox_2(t)*sqrt(2*(h_liq-h_2(t)));
        
%         m_dot_ox_in(t) = C_inj * (G_SPI(t)); % + G_HEM(t))/2;
%         
%         del_m_ox_tank = -m_dot_ox_in(t)*del_time;
%         del_U_tank = (-m_dot_ox_in(t)*h_liq + Q_dot_tank)*del_time;
        
        if x_tank(t) > 1
            burn_time = t*del_time; % Record burn time to be the time at which the tank runs out of oxidizer
%             liq_gas_transition = 1;
        end
    else
        % OXIDIZER TANK CALCULATIONS FOR VAPOUR ONLY %
        rho_tank(t) = m_ox_tank(t)./V_tank;
        u_tank(t) = U_tank(t)/m_ox_tank(t);
        uin.rho = rho_tank(t);
        uin.u = u_tank(t);
        if abs(uerror(T_tank(t), uin)) > u_tank_eps
            T_tank(t) = secant(@uerror, T_tank(t), uin);
        end
        
        thermo_span_tank = thermoSpanWagner(rho_tank(t), T_tank(t), {'p','h'});
        p_tank(t) = thermo_span_tank(1);
        h_tank(t) = thermo_span_tank(2) + 7.3397e+05; % Convert from Span-Wagner enthalpy convention to NIST
        h_discharge = h_tank(t);
        
%         if liq_gas_transition % 1-step transition from gas to liquid to smooth things over
%             rho_discharge = (6*rho_tank(t) + rho_liq)/7;
%             liq_gas_transition = 0;
%         else
            rho_discharge = rho_tank(t);
%         end
        
    end
    
%     if p_cc(t) > p_tank(t) % Check that the chamber pressure is sufficiently lower than the tank pressure
%         error('MuleSim2:highPressure', 'The combustion chamber pressure is too high (p_cc=%7.0f Pa, p_tank=%7.0f Pa, ) \n', p_cc(t), p_tank(t))
%     end
    
    % Iteratively solve for mass flow
%     m = 0;
%     m_dot_ox_temp = 0;
%     p_stag_temp = 0;
%     while (abs(m_dot_ox_temp - m_dot_ox_in(t)) > 0.1*m_dot_ox_in(t) || abs(p_stag_temp-p_stag(t)) > 0.1*p_stag(t)) && m < 100 % Until m_dot_ox_in and p_stag converges
        
        if (p_tank(t)-p_cc(t)-p_feed) < 0 % Check that the chamber pressure is sufficiently lower than the tank pressure
            error('MuleSim2:highPressure', 'The combustion chamber pressure is too high (p_cc=%7.0f Pa, p_tank=%7.0f Pa, ) \n', p_cc(t), p_tank(t))
        end
        
        % Might have different C_inj for a gas vs. a liquid
        m_dot_ox_in(t) = C_inj*sqrt(2*rho_discharge*(p_tank(t)-p_cc(t)-p_feed)); % Incompressible fluid assumption (better than nothing)
%         m_dot_ox_temp = m_dot_ox_in(t);
        
        % COMBUSTION CHAMBER CALCULATIONS %
        A_cc(t) = pi()*r_cc(t)^2; % Port geometry
        G(t) = m_dot_ox_in(t)/A_cc(t); % Mass flux
        del_r_cc(t) = a*G(t)^n * del_time; % Regression rate law
        m_dot_f(t) = 2*pi()*r_cc(t)*L * rho_f * (del_r_cc(t)/del_time); % Fuel mass flow rate from geometry

        % Iteratively solve for m_dot_f
        k = 0;
        m_dot_f_temp = 0;
        while abs(m_dot_f_temp - m_dot_f(t)) > 0.01*m_dot_f(t) && k < 100 % until m_dot_f converges
            m_dot_f_temp = m_dot_f(t);
            m_dot_cc(t) = m_dot_f(t) + m_dot_ox_in(t); % Total mass flow in
            G(t) = (m_dot_ox_in(t) + m_dot_cc(t))/(2*A_cc(t)); % Average mass flux accross entire port
            del_r_cc(t) = a*G(t)^n * del_time; % Regression rate law using updated G
            m_dot_f(t) = 2*pi()*r_cc(t)*L * rho_f * (del_r_cc(t)/del_time); % Fuel mass flow rate from geometry
            k = k+1;
        end
        OF(t) = m_dot_ox_in(t)./m_dot_f(t); % Calculate Oxidizer-Fuel Ratio

        % Iteratively solve for p_stag
%         k = 0;
%         p_stag_temp = 0;
%         while abs(p_stag_temp-p_stag(t)) > 0.01*p_stag(t) && k < 100 % Iterate until p_stag converges
            p_stag_temp = m_dot_cc(t)/(zeta_d*A_th)*sqrt(T_stag(t)*R_cc(t)/k_cc(t)*((k_cc(t)+1)/2)^((k_cc(t)+1)/(k_cc(t)-1))); % Steady state choked flow expression for pressure
            
            p_cc(t) = p_stag_temp*((T_cc(t)/T_stag(t)).^(k_cc(t)/(k_cc(t)-1))); % Velocity correction for stagnation pressure

            thermo_comb = CEAProp(OF(t), p_cc(t), {'T', 'rho', 'cp', 'k', 'M'});
            thermo_comb = num2cell(thermo_comb);
            [T_stag(t), rho_cc(t), cp_cc(t), k_cc(t), M_cc(t)] = thermo_comb{:};
            R_cc(t) = R_const/M_cc(t);

            T_stag(t) = T_stag(t)*zeta_cstar^2;

            vel_cc(t) = m_dot_cc(t)/(rho_cc(t)*A_cc(t));
            T_cc(t) = T_stag(t) - vel_cc(t)^2/(2*cp_cc(t));
            p_stag(t) = p_cc(t)*(T_stag(t)/T_cc(t)).^(k_cc(t)/(k_cc(t)-1)); % Re-calculate p_stag using new T_cc, T_stag, and k_cc
            
%             k = k + 1;
%         end
        
%         m_dot_ox_in(t) = C_inj*sqrt(2*rho_discharge*(p_tank(t)-p_cc(t)-p_feed)); % Incompressible fluid assumption
%         m = m+1;
%     end
    
    % NOZZLE AND THRUST CALCULATIONS %
    Ain.k = k_cc(t);
    Ain.A = A_ratio_nozzle;
    if abs(Aerror(Ma(t), Ain)) > A_ratio_nozzle_eps % Nozzle mach number solver
        Ma(t) = secant(@Aerror, Ma(t), Ain);
    end
    p_exit(t) = p_stag(t)./(1+(k_cc(t)-1)./2.*Ma(t).^2).^(k_cc(t)/(k_cc(t)-1)); % (Pa) Pressure of gas at exit of nozzle
    if x_tank(t) >= 1 % Physically very wrong, but the condition is used to hide the transition that results from overexpansion
        p_exit(t) = p_stag(t)*(2/(k_cc(t)+1))^(k_cc(t)/(k_cc(t)-1)); % Throat pressure
        vel_exit(t) = sqrt(2*k_cc(t)*R_cc(t)*T_stag(t)/(k_cc(t)+1)); % Throat velocity
        A_ratio_nozzle_eff = 1; % Throat area ratio
    else
        T_exit(t) = T_stag(t)./(1+(k_cc(t)-1)./2.*Ma(t).^2); % (K) Temperature of gas at exit of nozzle
        vel_exit(t) = Ma(t).*sqrt(k_cc(t).*R_cc(t).*T_exit(t)); % (m/s) Velocity of gas at exit of nozzle
        A_ratio_nozzle_eff = A_ratio_nozzle; % Effective nozzle area ratio is actual nozzle area ratio
    end
    F_thrust(t) = zeta_CF*(m_dot_cc(t).*vel_exit(t) + (p_exit(t)-p_atm).*A_th.*A_ratio_nozzle_eff); % (N) Rocket Motor Thrust!!!
    
    % ITERATE FORWARD IN TIME %
    if t < time_max/del_time && m_ox_tank(t)/m_ox_tank_init > 0.03
    % If less than max burn time and more than 3% of oxidizer is left, step forward in time
        % Oxidizer Tank Values
        del_m_ox_tank = -m_dot_ox_in(t)*del_time;
        del_U_tank = (-m_dot_ox_in(t)*h_discharge + Q_dot_tank)*del_time;
        m_ox_tank(t+1) = m_ox_tank(t) + del_m_ox_tank;
        U_tank(t+1) = U_tank(t) + del_U_tank;
        x_tank(t+1) = x_tank(t);
        T_tank(t+1) = T_tank(t);
        m_dot_ox_in(t+1) = m_dot_ox_in(t);
        % Combustion Chamber Values
        m_tot(t+1) = m_tot(t) - m_dot_cc(t)*del_time;
        r_cc(t+1) = r_cc(t) + del_r_cc(t);
        p_cc(t+1) = p_cc(t);
        p_stag(t+1) = p_stag(t);
        T_cc(t+1) = T_cc(t);
        T_stag(t+1) = T_stag(t);
        R_cc(t+1) = R_cc(t);
        k_cc(t+1) = k_cc(t);
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
%%%% CHANGE CALCULATIONS TO ONLY USE THE LIQUID BURN PORTION %%%%
r_dot = del_r_cc./del_time; % (m/s)
d_port_f = 2*r_cc(end); % (m)
m_out = trapz(m_dot_cc)*del_time; % (kg)
m_f = trapz(m_dot_f)*del_time; % (kg)
m_ox = m_out-m_f; % (kg)
I_tot = trapz(F_thrust)*del_time; % (N*s)
v_e = I_tot/m_out; % (m/s)
I_sp = v_e/9.81; % (s)
c_star = mean(p_cc)*A_th/mean(m_dot_cc); % (m/s)


% Save metadata of simulation
time_sim = toc; % Save how long the simulation took (minus loading data and plotting figures)
time_stamp = datestr(datetime); % Save the current time and date
computer_name = computer; % Save which computer the code is running on
function_name = mfilename; % Save the name of the function that is being called

TCurve = [I_tot, mean(F_thrust), m_out];
TCurve = vertcat(TCurve, [time; F_thrust; m_tot]'); % Output thrust and mass vs time for trajectory analysis


%% PLOTS %%
if graph == 1 % Graph if indicated
    
%     X_actual = valid.Motor_Validation.UofTDeliveranceIITankPressure;
%     X_actual = rmmissing(X_actual); % Remove trailing NaNs
%     time_actual = (0:0.01:(length(X_actual)-1)*0.01)';
%     
%     figure
%     hold on
%     plot(time, p_tank./1e6)
%     plot(time_actual, X_actual./1e6)
%     title('University of Toronto - Boundless')
%     xlabel('Time (s)')
%     ylabel('Tank Pressure (MPa)')
%     legend('MuleSim3', 'Test Data')
% %     axis([0 20 0 6000])
% %     saveas(gcf, ['Output\DeliveranceII_ptank_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off
    
    X_actual = valid.Motor_Validation.Pheonix1AStaticChamberPressure;
    X_actual = rmmissing(X_actual); % Remove trailing NaNs
    time_actual = (0:0.01:(length(X_actual)-1)*0.01)';
%     X2_actual = valid.Motor_Validation.Pheonix1AModelChamberPressure;
%     X2_actual = rmmissing(X2_actual); % Remove trailing NaNs
%     time_actual_2 = (0:0.01:(length(X2_actual)-1)*0.01)';
    
    figure
    hold on
    plot(time, p_cc./1e6)
%     plot(time_actual_2, X2_actual./1e6)
%     plot(time_actual, X_actual./1e6)  
    plot(time, p_tank./1e6)
    title('University of Kwazulul-Natal - Phoenix 1A')
    xlabel('Time (s)')
    ylabel('Chamber Pressure (MPa)')
    legend('MuleSim3', 'Test Data')
%     axis([0 20 0 6000])
%     saveas(gcf, ['Output\Phoenix1A_pcc_', datestr(now, 'yyyy-mm-dd'), '.png']);
    hold off
    
    X_actual = valid.Motor_Validation.Pheonix1AStaticThrust;
    X_actual = rmmissing(X_actual); % Remove trailing NaNs
    time_actual_1 = (0:0.01:(length(X_actual)-1)*0.01)';
%     X2_actual = valid.Motor_Validation.Phase2StanfordThrust;
%     X2_actual = rmmissing(X2_actual); % Remove trailing NaNs
%     time_actual_2 = (0:0.01:(length(X2_actual)-1)*0.01)';
    
    figure
    hold on
    plot(time, F_thrust)
    plot(time_actual_1, X_actual)
%     plot(time_actual_2, X2_actual)
    title('University of Kwazulul-Natal - Phoenix 1A')
    xlabel('Time (s)')
    ylabel('Thrust (N)')
    legend('MuleSim3', 'Test Data')
%     axis([0 20 0 6000])
%     saveas(gcf, ['Output\Phoenix1A_Fthrust_', datestr(now, 'yyyy-mm-dd'), '.png']);
    hold off

end


%% CALCULATE ERROR %%
% if exist('F_act') % if F_act
%     MAE = mean(abs(F_act'-F_thrust_error));
%     e = MAE/mean(F_act)*100;
% end
    

%% SAVE OUTPUT TO SPREADSHEET %%
if save == 1 % Save if indicated
    [file,path] = uiputfile(['Output\', mfilename, '.xlsx']); % User input file name
    if file ~= 0
        filename = fullfile(path, file);
        % Record fuel and oxidizer data
        ox_k = fieldnames(CEA.OX); % Read oxidizer info from CEA file
        for k = 1:length(ox_k)
            Parameter_OX{k} = ['Oxidizer ', num2str(k), ' Mass Fraction'];
            Symbol_OX{k} = ox_k{k};
            Value_OX{k} = CEA.OX.(ox_k{k}).wtfrac;
            Units_OX{k} = '-';
        end
        fuel_k = fieldnames(CEA.FUEL); % Read fuel info from CEA file
        for k = 1:length(fuel_k)
            Parameter_FUEL{k} = ['Fuel ', num2str(k), ' Mass Fraction'];
            Symbol_FUEL{k} = fuel_k{k};
            Value_FUEL{k} = CEA.FUEL.(fuel_k{k}).wtfrac;
            Units_FUEL{k} = '-';
        end

        % Add output data to table
        Parameter_out = {'Burn Time'; 'Peak Thrust'; 'Average Thrust'; ...
            'Total Impulse'; 'Specific Impulse'; 'Average O/F Ratio'; ...
            'Average Oxidizer Mass Flow'; 'Average Fuel Mass Flow'; ...
            'Total Oxidizer Mass Consumed'; 'Total Fuel Mass Consumed'; ...
            'Average Regression Rate'; 'Final Port Diameter'; ...
            'Maximum Combustion Chamber Pressure'; ...
            'Maximum Combustion Chamber Temperature'}; 
        Symbol_out = {'burn_time'; 'max(F_thrust)'; 'mean(F_thrust)'; 'I_tot'; ...
            'I_sp'; 'mean(OF)'; 'mean(m_dot_ox_in)'; 'mean(m_dot_f)'; 'm_ox'; ...
            'm_f'; 'mean(r_dot)'; 'd_port_f'; 'max(T_cc)'; 'max(p_cc)'};
        Value_out = {burn_time; max(F_thrust); mean(F_thrust); I_tot; I_sp; ...
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
        warning('off','MATLAB:xlswrite:AddSheet'); % Stop MATLAB from complaining when adding a new sheet
        writetable(output, filename, 'Sheet', 'Summary');

        % Output Time-series data
        output_time = table(['[s]'; num2cell(time')], ['[N]'; num2cell(F_thrust')], ...
            ['[K]'; num2cell(T_cc')], ['[Pa]'; num2cell(p_cc')], ['[-]'; num2cell(OF')], ...
            ['[K]'; num2cell(T_tank')], ['[Pa]'; num2cell(p_tank')], ...
            ['[-]'; num2cell(x_tank')]);
        output_time.Properties.VariableNames = {'Time', 'Thrust', ...
            'T_cc', 'p_cc', 'OF', 'T_tank', 'p_tank', 'x_tank'};
        writetable(output_time, filename, 'Sheet', 'Time Data');
        
        % Delete extra sheets that get created with a new excel file
        objExcel = actxserver('Excel.Application');
        objExcel.Workbooks.Open(filename);
        try  % Throws an error if the sheets do not exist
              objExcel.ActiveWorkbook.Worksheets.Item(['Sheet1']).Delete; % Delete sheets
              objExcel.ActiveWorkbook.Worksheets.Item(['Sheet2']).Delete;
              objExcel.ActiveWorkbook.Worksheets.Item(['Sheet3']).Delete;
        catch
              % Do nothing
        end
        objExcel.ActiveWorkbook.Save; % Save, close, and clean up
        objExcel.ActiveWorkbook.Close;
        objExcel.Quit;
        objExcel.delete;
    end
end

%% SUBFUNCTIONS %%
    function [out] = CEAProp(OF, p, in)
    %MASSFRAC Calculates the mass fractions of products for a given reaction
    %   'OF' is the oxidizer to fuel weight ratio
    %   'p' is the pressure of the combustion chamber
    %   'alpha' is an array of weight fractions
        if OF < CEA.OF(1) % Check that query is inside bounds
            error('CEAProp:lowOF', 'Outside O/F range: \n OF = %f', OF)
        elseif OF > CEA.OF(end)
            error('CEAProp:highOF', 'Outside O/F range: \n OF = %f', OF)
        elseif p < CEA.p(1)
            error('CEAProp:lowP', 'Outside pressure range: \n p = %f Pa', p)
        elseif p > CEA.p(end)
            error('CEAProp:highP', 'Outside pressure range: \n p = %f Pa', p)
        end
        for mm = 1:length(CEA.OF)
            if OF < CEA.OF(mm) % Find first index that is larger than OF
                break
            end
        end
        for nn = 1:length(CEA.p)
            if p < CEA.p(nn) % Find first index that is larger than p
                break
            end
        end
        out = zeros(size(in));
        for kk = 1:length(in) % Iterate through all products
            % Do some bilinear interpolation (equations from wikipedia)
            OF_int = [CEA.OF(mm)-OF, OF-CEA.OF(mm-1)];
            p_int = [CEA.p(nn)-p; p-CEA.p(nn-1)];
            C_int = 1/((CEA.OF(mm)-CEA.OF(mm-1))*(CEA.p(nn)-CEA.p(nn-1)));
            out_int = CEA.(in{kk})(mm-1:mm,nn-1:nn);
            out(kk) = C_int.* OF_int*out_int*p_int;
        end
    end


    function [val_out] = thermoSat(val_in, in, out)
    %THERMOSAT returns thermodynamic properties at saturation
    %   'val_in' is a double specifying the value of thermodynamic property inputted
    %   'in' is a string specifying the given input thermodynamic property 
    %   'out' as a cell array specifying the desired output thermodynamic property
    %   'val_out' is a double specifying the value of thermodynamic property returned
        val_out = zeros(size(out));
        for mm = 1:length(out) % For each desired output
            % Locate which variable column is input, and which is output
            ii = find( (in == N2Osat.meta(1,:)), 1, 'first');
            jj = find( (out(mm) == N2Osat.meta(1,:)), 1, 'first');
            % Interpolate between the two values
            for kk = 1:length(N2Osat.data(:,ii))
                if val_in < N2Osat.data(kk,ii)
                    break
                end
            end
            val_out(mm) = (val_in - N2Osat.data(kk-1,ii))./(N2Osat.data(kk,ii)-N2Osat.data(kk-1,ii)).*(N2Osat.data(kk,jj)-N2Osat.data(kk-1,jj)) + N2Osat.data(kk-1,jj);
        end
    end


    function [out] = thermoSpanWagner(rho, T, in)
    %THERMOSPANWAGNER Calculates thermodynamic properties for N2O as a non-ideal gas
    %   'in' is a cell array containing the desired output parameters
    %   'rho' is the density in kg/m^3
    %   'T' is the temperature in K
    %   'out' is an array containing the output values in the order listed in 'in'

        % Hardcode in data for N2O (from "Modeling Feed System Flow Physics for Self-Pressurizing Propellants)
        R = 8.3144598/44.0128*1000; % (J/kg*K) Gas constant
        T_c = 309.52; % (K) Critical Temperature
        rho_c = 452.0115; % (kg/m^3) Critical Density 
        n0 = [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, ...
            -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038];
        n1 = n0(1:5); n2 = n0(6:12);
        a1 = 10.7927224829;
        a2 = -8.2418318753;
        c0 = 3.5;
        v0 = [2.1769, 1.6145, 0.48393];
        u0 = [879, 2372, 5447];
        t0 = [0.25, 1.125, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5];
        d0 = [1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2];
        P0 = [1, 1, 1, 2, 2, 2, 3];
        t1 = t0(1:5); t2 = t0(6:12);
        d1 = d0(1:5); d2 = d0(6:12);

        % Calculate non-dimensional variables
        tau = T_c/T;
        delta = rho/rho_c;

        % Calculate explicit helmholtz energy and derivatives (from 
        ao = a1 + a2*tau + log(delta) + (c0-1)*log(tau) + sum(v0.*log(1-exp(-u0.*tau./T_c)));
        ar = sum(n1.*tau.^t1.*delta.^d1) + sum(n2.*tau.^t2.*delta.^d2.*exp(-delta.^P0));
        ao_tau = a2 + (c0-1)/tau + sum(v0.*u0./T_c.*exp(-u0.*tau./T_c)./(1-exp(-u0.*tau./T_c)));
        ao_tautau = -(c0-1)/tau.^2 + sum(-v0.*u0.^2./T_c.^2.*exp(-u0.*tau./T_c)./(1-exp(-u0.*tau./T_c)).^2);
        ar_tau = sum(n1.*t1.*tau.^(t1-1).*delta.^d1) + sum(n2.*t2.*tau.^(t2-1).*delta.^d2.*exp(-delta.^P0));
        ar_tautau = sum(n1.*t1.*(t1-1).*tau.^(t1-2).*delta.^d1) + sum(n2.*t2.*(t2-2).*tau.^(t2-2).*delta.^d2.*exp(-delta.^P0));
        ar_delta = sum(n1.*d1.*delta.^(d1-1).*tau.^t1) + sum(n2.*tau.^t2.*delta.^(d2-1).*(d2-P0.*delta.^P0).*exp(-delta.^P0));
        ar_deltadelta = sum(n1.*d1.*(d1-1).*delta.^(d1-2).*tau.^t1) + sum(n2.*tau.^t2.*delta.^(d2-2).*((d2-P0.*delta.^P0).*(d2-1-P0.*delta.^P0)-P0.^2.*delta.^P0).*exp(-delta.^P0));
        ar_deltatau = sum(n1.*d1.*t1.*delta.^(d1-1).*tau.^(t1-1)) + sum(n2.*t2.*tau.^(t2-1).*delta.^(d2-1).*(d2-P0.*delta.^P0).*exp(-delta.^P0));

        out = zeros(size(in));
        for kk = 1:length(in)
            switch in{kk}
                case 'p' % Pressure (Pa)
                    out(kk) = rho*R*T*(1+delta*ar_delta);
                case 'u' % Specific internal energy (J/kg)
                    out(kk) = R*T*tau*(ao_tau+ar_tau);
                case 's' % Specific entropy (J/kg*K)
                    out(kk) = R*(tau*(ao_tau+ar_tau)-ao-ar);
                case 'h' % Specific enthalpy (J/kg)
                    out(kk) = R*T*(1+tau*(ao_tau+ar_tau)+delta*ar_delta);
                case 'cv' % Specific heat constant pressure (J/kg*K)
                    out(kk) = R*-tau^2*(ao_tautau+ar_tautau);
                case 'cp' % Specific heat constant pressure (J/kg*K)
                    out(kk) = R*(-tau^2*(ao_tautau+ar_tautau) + (1+delta*ar_delta-delta*tau*ar_deltatau)^2/(1+2*delta*ar_delta+delta^2*ar_deltadelta));
                case 'a' % Speed of sound (m/s)
                    out(kk) = sqrt(R*T*(1+2*delta*ar_delta+delta^2*ar_deltadelta - (1+delta*ar_delta-delta*tau*ar_deltatau)^2/(tau^2*(ao_tautau+ar_tautau))));
                otherwise
                    error('Invalid input')
            end
        end
    end

%     function [y_np1] = rungeKutta4(
%         
%     end


    function [x] = secant(fun, x1, in)
    %SECANT is a zero-finding function based on the secant method
    %   'fun' is the function handle for which the zero is desired
    %   'x1' is the initial guess
    %   'in' is a struct containing any additional inputs 'fun' might require
    %   'x' is the value of the zero
        x_eps = x1*0.005; % Set the tolerance to be 0.5% of initial guess
        x2 = x1-x1*0.01; % Set a second point 1% away from the original guess
        F1 = fun(x1, in); % Evaluate function at x1
        F2 = fun(x2, in); % Evaluate function at x2
        kk = 1; % Set up counter
        kk_max = 1000;
        while abs(x2-x1)>=x_eps && kk<kk_max % While error is too large and counter is less than max
            x3 = x2 - F2*(x2-x1)/(F2-F1);
            x1 = x2; % Move everything forward
            x2 = x3;
            F1 = F2;
            F2 = fun(x2, in);
            kk = kk+1;
        end
        x = x2;
    end


    % SET UP FUNCTIONS TO BE ITERATIVELY SOLVED %
    function V = Verror(T, in) % Finds the difference between the estimated and actual tank volume
        thermo = thermoSat(T, 'T', {'rho_liq','rho_vap','u_liq','u_vap'});
        thermo = num2cell(thermo);
        [rho_l,rho_v,u_l,u_v] = thermo{:};
        x = (in.U/in.m - u_l)/(u_v - u_l);
        V = in.m*((1-x)/rho_l + x/rho_v) - in.V;
    end
    
    function U = uerror(T, in) % Finds the difference between the estimated and actual tank internal energy
        U = (thermoSpanWagner(in.rho, T, {'u'}) + 7.3397e+5) - in.u;
    end

    function A = Aerror(M, in) % Finds the difference between the estimated and actual nozzle ratio
    	A = (1/M^2)*(2./(in.k+1).*(1+(in.k-1)./2.*M.^2)).^((in.k+1)./(in.k-1)) - in.A^2;
    end

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

