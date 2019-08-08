function [p] = SimShell(CEA_f, input_f)
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
    CEA_f = 'MuleSim3CEA';
elseif nargin == 1
    input_f = 'MuleSim3INPUT';
elseif nargin > 2
    error('Invalid number of inputs');
end
input = readtable([input_f, '.xlsx']);
CEA = load([CEA_f, '.mat']);
N2Osat = load('N2Osat.mat');

% Pre-declare variables before they're poofed into existence by eval
[del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, d_f, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    zeta_d, zeta_cstar, zeta_CF, d_th, p_atm, A_ratio_nozzle, graph, save] = deal(NaN);

for k = 1:length(input.Symbol)
    if ~exist(input.Symbol{k}, 'var')
        error('Variable ''%s'' in the input spreadsheet is not declared in the MATLAB code', input.Symbol{k})
    end
    eval([input.Symbol{k} '= input.Value(' num2str(k) ');']); % Read input data from spreadsheet
end

design = [del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, d_f, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    zeta_d, zeta_cstar, zeta_CF, d_th, p_atm, A_ratio_nozzle];

p.design = design;
p.CEA = CEA;
p.N2Osat = N2Osat;

%% SIMULATION

p = Surrocket(p);

%% PLOTS %%
if graph == 1 % Graph if indicated
%     figure
%     hold on
%     plot(time, p_cc./1e+05) % (bar) 6894.76 Pa per psi
%     title('Pheonix 1A Hot-Fire Test Chamber Pressure')
%     xlabel('Time (s)')
%     ylabel('Pressure (bar)')
%     axis([0 18 0 40])
% %     saveas(gcf, ['Output\Pheonix_Pcc_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off
     
%     figure
%     hold on
%     plot(time, F_thrust) % (N)
%     title('Pheonix 1A Hot-Fire Test Thrust')
%     xlabel('Time (s)')
%     ylabel('Thrust (N)')
%     axis([0 18 0 4000])
% %     saveas(gcf, ['Output\Pheonix_Thrust_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off

%     figure
%     hold on
%     plot(time, F_thrust) % (N)
%     title('3" Diameter Flight Data Thrust')
%     xlabel('Time (s)')
%     ylabel('Thrust (N)')
%     axis([0 10 -200 2000])
% %     saveas(gcf, ['Output\3inch_Stanford_Thrust_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off

%     figure
%     hold on
%     plot(time, F_thrust) % (N)
%     title('Phase 2 Simulation Comparison')
%     xlabel('Time (s)')
%     ylabel('Thrust (N)')
%     axis([0 16 0 6000])
% %     saveas(gcf, ['Output\Phase_2_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off
    
%     figure
%     hold on
%     plot(time, F_thrust/9.80665) % (kgf)
%     title('RATTWORKS M-900 Thrust')
%     xlabel('Time (s)')
%     ylabel('Thrust (kgf)')
%     axis([0 10.5 0 112])
% %     saveas(gcf, ['Output\RATTWORKS_M-900_Thrust_', datestr(now, 'yyyy-mm-dd'), '.png']);
%     hold off

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
%     saveas(gcf, ['Output\Boundless_', datestr(now, 'yyyy-mm-dd'), '.png']);
    hold off
end

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

