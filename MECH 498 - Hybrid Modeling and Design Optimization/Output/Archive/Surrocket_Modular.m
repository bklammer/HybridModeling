function [out, p] = Surrocket(p)
%SURROCKET Models a hybrid rocket from design to apogee
%   Detailed explanation goes here


%% INPUTS %%
% Load input files
if nargin == 0 % Set default files to read from in case no input
    input_f = 'MuleSim4INPUT';
    CEA_f = 'MuleSim3CEA';
elseif nargin == 1
    input_f = 'MuleSim4INPUT';
elseif nargin > 2
    error('Invalid number of inputs');
end
input = readtable([input_f, '.xlsx']);
CEA = load([CEA_f, '.mat']);
N2Osat = load('N2Osat.mat');

% Pre-declare variables before they're poofed into existence by eval
[del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, m_f_init, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
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
[p.T_tank, rho_liq, rho_vap, u_liq, u_vap] = thermo_sat_init{:};
p.x_tank = (V_tank/m_ox_tank_init - 1/rho_liq)/(1/rho_vap - 1/rho_liq); % Calculate vapour mass fraction
u_tank = (p.x_tank*u_vap + (1-p.x_tank)*u_liq); % (J/kg*K) Calculate specific internal energy
p.U_tank = m_ox_tank_init*u_tank; % (J) Calculate total internal energy

% Calculation of constants
p.C_inj = n_inj * Cd_inj * (pi()*(d_inj/2)^2); % (m^2) Injector coefficient
p.A_th = pi()*(d_th/2)^2; % (m^2) Nozzle Throat Area
p.V_tank_eps = 0.001*V_tank; % (m^3) Set acceptable tank volume error to 0.1% of tank volume
p.u_tank_eps = 0.001*u_tank; % (m^3) Set acceptable tank specific internal energy error to 0.1% of initial tank specific internal energy
p.A_ratio_nozzle_eps = 0.001*A_ratio_nozzle; % (m^3) Set acceptable nozzle area ratio error to 0.1% of nozzle area ratio

% Setting initial conditions
p.m_ox_tank = m_ox_tank_init; % (kg) Mass of oxidizer in tank
p.r_cc = d_port_init/2; % (m) Combustion chamber port radius
p.p_cc = p_cc_init; % (kg) Array storing mass of each product
p.T_cc = T_cc_init; % (K) Combustion chamber temperature
p.T_stag = T_cc_init;
p.R_cc = 8.3144598/0.02897; % (J/kg*K) Set initial fluid in chamber to be air
p.k_cc = 1.4; % Values for air at 298K
p.Ma = 3; % Initialize the mach number to three

p.V_tank = V_tank; 
p.del_time = del_time;
p.t_max = time_max/del_time;
p.N2Osat = N2Osat;
p.CEA = CEA;
p.a = a;
p.n = n;
p.zeta_d = zeta_d;
p.zeta_cstar = zeta_cstar;
p.A_ratio_nozzle = A_ratio_nozzle;
p.zeta_CF = zeta_CF;
p.p_atm = p_atm;
p.rho_f = rho_f;
p.L = L;


%% MODEL CALCS

p.t = 1;
% Hybrid Motor
while p.m_ox_tank(p.t)/p.m_ox_tank(1) > 0.03 && p.t < p.t_max % Hybrid Rocket Model
    p = SelfPressurizedOxidizerTank(p);
    p = CombustionChamber(p);
    p = Nozzle(p);
    p.t = p.t + 1;
end
time_sim = toc
time = 0:del_time:(del_time*(p.t-1)); % Create a time vector
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
time_sim = toc

% Structural Estimation
% 
% 
% % Trajectory Analysis
% x.t = 0;
% while x.vel_rocket > 0 && t < t_max
%     x = trajectory(x);
% end
% 
% out = x.apogee;
% 

if graph == 1
    figure(1)
    hold on
    plot(time, p.F_thrust./4.44822)
    plot(time, p.p_tank./6894.76)
    plot(time, p.p_cc./6894.76)
    title('Boundless - University of Washington Results')
    xlabel('Time (s)')
    ylabel('Thrust (lbf) and Pressure (psi)')
    legend('Thrust', 'Tank Pressure', 'Chamber Pressure')
    axis([-1 20 -50 1400])
%     saveas(gcf, ['Output\Boundless_', datestr(now, 'yyyy-mm-dd'), '.png']);
    hold off
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

end




                        
% 
%                                         @@&              
%                                  %%   @  @@  @           
%                              @      *@     @             
%                    @@@@%   @      @@     @/#@@           
%                  .@@@     @ /@  @,            @          
%                  .@@@     (   @              %@          
%                   @@@@    @       @@   @    %@           
%                     @@@@ .%,@.  @ *   @    @             
%                       @@@@( @ @@ @   /  @/   @@          
%                      @  *@@@@@     #         @@/         
%                     @       @@@@@@@@.     ,@@@@          
%                     @     @      %@@@@@@@@@@@.           
%                    & @(                                  
% 



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
