function [e] = OxTankValidation(input_f)
%TANKSIM Model of self-pressurized propellant tank
% 'input_f' is a string containing the name of a .xlsx file containing tank data
% All input files must be in the same folder as the 'MuleSim3.m' or have the folder path included

% Benjamin Klammer - 2019
% All calculations performed in metric units

%% INPUTS %%
% Load input files
if nargin == 0 % Set default files to read from in case no input
    input_f = 'TankINPUT';
elseif nargin > 1
    error('Invalid number of inputs');
end
input = readtable([input_f, '.xlsx']);
N2Osat = load('N2Osat.mat');
p = load('p_tank_validation.mat');

tank_names = input.Properties.VariableNames;
for k = 4:length(tank_names) % Simulate all the different tanks back-to-back
tic; % Start recording the simulation runtime

tank_input = num2cell(input.(tank_names{k}));
[V_tank, m_ox_tank_init, p_tank_init, C_inj, p_atm, del_time, time_max] = tank_input{:};

%% INITIAL CALCULATIONS %%
% Calculate initial thermodynamic properties from tank pressure
thermo_sat_init = thermoSat(p_tank_init, 'p', {'T', 'rho_liq', 'rho_vap', 'u_liq', 'u_vap'});
thermo_sat_init = num2cell(thermo_sat_init);
[T_tank, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_init{:};
x_tank = (V_tank/m_ox_tank_init - 1/rho_liq)/(1/rho_vap - 1/rho_liq); % Calculate vapour mass fraction
u_tank = (x_tank*u_vap + (1-x_tank)*u_liq); % (J/kg*K) Calculate specific internal energy
U_tank = m_ox_tank_init*u_tank; % (J) Calculate total internal energy

% Calculation of constants and initial conditions
Q_dot_tank = 0; % (W) Ignore the effect of heat losses to the tank walls (for now)
V_tank_eps = 0.001*V_tank; % (m^3) Set acceptable tank volume error to 0.1% of tank volume
u_tank_eps = 0.001*u_tank; % (m^3) Set acceptable tank specific internal energy error to 0.1% of initial tank specific internal energy
m_ox_tank = m_ox_tank_init; % (kg) Mass of oxidizer in tank

% Initialize all other variables so that they don't get overwritten by
% previous longer burns
p_tank = 0;


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
        
        % Thermo properties downstream of injector for HEM model
        s_liq = thermoSat(T_tank(t), 'T', {'s_liq'});
        thermo_sat_inj = thermoSat(p_atm, 'p', {'s_liq', 's_vap', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap'});
        thermo_sat_inj = num2cell(thermo_sat_inj);
        [s_liq_2, s_vap_2, h_liq_2, h_vap_2, rho_liq_2, rho_vap_2] = thermo_sat_inj{:};
        x_2(t) = (s_liq - s_liq_2)/(s_vap_2 - s_liq_2); % Mass fraction downstream assuming isentropic expansion at saturation (only use for HEM mass flow)
        rho_ox_2(t) = 1/(x_2(t)/rho_vap_2 + (1 - x_2(t))/rho_liq_2);
        h_2(t) = x_2(t)*h_vap_2 + (1 - x_2(t))*h_liq_2;

        % Calculate mass flow through the injector
        G_SPI(t) = sqrt(2*rho_liq*(p_tank(t)-p_atm));
        G_HEM(t) = rho_ox_2(t)*sqrt(2*(h_liq-h_2(t)));
        m_dot_ox_in(t) = C_inj * (G_SPI(t) + G_HEM(t))/2;
        
        del_m_ox_tank = -m_dot_ox_in(t)*del_time;
        del_U_tank = (-m_dot_ox_in(t)*h_liq + Q_dot_tank)*del_time;
        
        if x_tank(t) > 1
            burn_time(k) = t*del_time; % Record burn time to be the time at which the tank runs out of oxidizer
            liq_gas_transition = 1;
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
        
        if liq_gas_transition % 1-step transition from gas to liquid to smooth things over
            rho_tank(t) = (6*rho_tank(t) + rho_liq)/7;
            liq_gas_transition = 0;
        end
        
        m_dot_ox_in(t) = C_inj*sqrt(2*rho_tank(t)*(p_tank(t)-p_atm)); % Incompressible fluid assumption (better than nothing)
        
        del_m_ox_tank = -m_dot_ox_in(t)*del_time;
        del_U_tank = (-m_dot_ox_in(t)*h_tank(t) + Q_dot_tank)*del_time;
    end
    
    % ITERATE FORWARD IN TIME %
    if t < time_max/del_time && m_ox_tank(t)/m_ox_tank_init > 0.03
    % If less than max burn time and more than 3% of oxidizer is left, step forward in time
        % Oxidizer Tank Values
        m_ox_tank(t+1) = m_ox_tank(t) + del_m_ox_tank;
        U_tank(t+1) = U_tank(t) + del_U_tank;
        x_tank(t+1) = x_tank(t);
        T_tank(t+1) = T_tank(t);
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

%% PLOTS %%

p_tank_actual = p.p_tank_validation.(tank_names{k});
p_tank_actual = rmmissing(p_tank_actual); % Remove trailing NaNs
time_actual = (0:0.01:(length(p_tank_actual)-1)*0.01)';

% figure
% hold on
% plot(time, p_tank./1e+06) % (MPa)
% plot(time_actual, p_tank_actual./1e+06) % (MPa)
% title(tank_names{k})
% xlabel('Time (s)')
% ylabel('Pressure (MPa)')
% axis([0 inf 0 5.5])
% legend('Model', 'Experiment')
% saveas(gcf, ['Output\' tank_names{k}, '_', datestr(now, 'yyyy-mm-dd'), '.png']);
% hold off

%% ERROR CALCULATION %%
index_LRO = burn_time(k)/del_time; % Only calculate error up to liquid-run-out
MAE = mean(abs(p_tank_actual(1:index_LRO)'-p_tank(1:index_LRO)));
e(k) = MAE/mean(p_tank_actual(1:index_LRO))*100;

%% OUTPUTS %%
p_tank_LRO(k) = p_tank(index_LRO);
m_dot_ox_avg(k) = mean(m_dot_ox_in(1:index_LRO));
time_sim(k) = toc; % Save how long the simulation took (minus loading data and plotting figures)
time_stamp{k} = datestr(datetime); % Save the current time and date
computer_name{k} = computer; % Save which computer the code is running on
function_name{k} = mfilename; % Save the name of the function that is being called


end % end for loop

%% SAVE OUTPUT TO SPREADSHEET %%

[file,path] = uiputfile(['Output\', mfilename, '.xlsx']); % User input file name
if file ~= 0
    filename = fullfile(path, file);

    % Add output data to table
    Parameter_out = {'Liquid Run Out Time'; 'Average Oxidizer Mass Flow'; ...
        'Final Tank Pressure'; 'Pressure Error'};
    Symbol_out = {'burn_time'; 'm_dot_ox_avg'; 'p_tank_LRO'; 'e'};
    Units_out = {'s'; 'kg/s'; 'Pa'; '%'};

    % Add metadata to table
    Parameter_meta = {'MATLAB Function Name'; 'Date and Time Simulated'; ...
        'Elapsed Simulation Time'; 'Operating System and Version'};
    Symbol_meta = {'function_name'; 'time_stamp'; 'time_sim'; 'computer_name'};
    Units_meta = {'-'; '-'; 's'; '-'};
    
    % Add values to table
    for k = 4:length(tank_names)
        Value_meta = {function_name{k}; time_stamp{k}; time_sim(k); computer_name{k}};
        Value_out = {burn_time(k); m_dot_ox_avg(k); p_tank_LRO(k); e(k)};
        Value_input = num2cell(input.(tank_names{k}));
        Value(:,k-3) = [Value_meta; Value_input; Value_out];
    end
    
    % Output to table
    Parameter = [Parameter_meta; input.Parameter; Parameter_out];
    Symbol = [Symbol_meta; input.Symbol; Symbol_out];
    Units = [Units_meta; input.Units; Units_out];
    PrinceFlight = Value(:,1);
    PrinceGround = Value(:,2);
    ZimmermanLow = Value(:,3);
    ZimmermanHigh = Value(:,4);
    VanPelt = Value(:,5);
    Zilliac = Value(:,6);
    
    % Combine data to table and write to file
    output = table(Parameter, Symbol, Units, PrinceFlight, PrinceGround, ...
        ZimmermanLow, ZimmermanHigh, VanPelt, Zilliac);
    warning('off','MATLAB:xlswrite:AddSheet'); % Stop MATLAB from complaining when adding a new sheet
    writetable(output, filename, 'Sheet', 'Summary');

%     % Output Time-series data
%     output_time = table(['[s]'; num2cell(time')], ['[Pa]'; num2cell(p_tank')]);
%     output_time.Properties.VariableNames = {'Time', 'p_tank'};
%     writetable(output_time, filename, 'Sheet', 'Time Data');

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


%% SUBFUNCTIONS %%
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

end
