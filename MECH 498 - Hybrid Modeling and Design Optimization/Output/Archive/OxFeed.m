function [] = OxFeed()
%OXFEED Prototype oxidizer feed system model

%% INPUTS

N2Osat = load('N2Osat.mat');
N2Osat = N2Osat.N2Osat;

t_max = 100000;
t = 1;
del_time = 0.001; % (s) time step

tank_eps = 0.000005; % (m^3) tank check volume must be this close to actual volume to progress
del_T_tank = 0.01; % (K) degrees of temperature to change on each iteration

V_tank = 0.012; % (m^3) From UofT, 12L tank
m_tank = 8; % (kg) From UofT
p_tank = 5648000; % (Pa) saturated N2O @ 298K

Cd_inj = 0.66; % from Dyer et al
r_inj = 0.0015; % (m) radius of injector
n_inj = 6; % number of injector orifices

p_cc = 3000000; % (Pa) combustion chamber pressure

Q_dot_in = 0; % Revise later on to inculde heat transfer

%% INITIAL CALCULATIONS
% Calculate initial thermodynamic properties from tank pressure
T_tank = thermoSat(p_tank, 'p', 'T'); 
h_liq = thermoSat(p_tank, 'p', 'h_liq');
h_vap = thermoSat(p_tank, 'p', 'h_vap');
rho_liq = thermoSat(p_tank, 'p', 'rho_liq');
rho_vap = thermoSat(p_tank, 'p', 'rho_vap');
u_liq = h_liq - p_tank/rho_liq;
u_vap = h_vap - p_tank/rho_vap;

x_tank = (V_tank/m_tank - 1/rho_liq)/(1/rho_vap - 1/rho_liq);
U_tank = m_tank*(x_tank*u_vap + (1-x_tank)*u_liq); % Calculate total internal energy

C_inj = n_inj * Cd_inj * (pi()*r_inj^2); % Injector coefficient


%% MODEL CALCULATIONS

while x_tank(t)<1 && t<t_max % while liquid is left in the chamber
    
    % Iteratively estimate tank temperature
    k = 0;
    while true
        x_tank(t) = (U_tank(t)/m_tank(t) - u_liq)/(u_vap - u_liq);
        V_tank_check = m_tank(t)*((1-x_tank(t))/rho_liq + x_tank(t)/rho_vap);
        
        if abs(V_tank-V_tank_check) < tank_eps
            break
        end
        
        if V_tank_check > V_tank
            T_tank(t) = T_tank(t) + del_T_tank;
        else
            T_tank(t) = T_tank(t) - del_T_tank;
        end
        
        % Recalculate thermo properties at new temperature T
        p_tank(t) = thermoSat(T_tank(t), 'T', 'p');
        h_liq = thermoSat(T_tank(t), 'T', 'h_liq');
        h_vap = thermoSat(T_tank(t), 'T', 'h_vap');
        rho_liq = thermoSat(T_tank(t), 'T', 'rho_liq');
        rho_vap = thermoSat(T_tank(t), 'T', 'rho_vap');
        u_liq = h_liq - p_tank(t)/rho_liq;
        u_vap = h_vap - p_tank(t)/rho_vap;
        
        k = k + 1;
        if k > 1000
            error('Tank volume did not converge in 1000 iterations');
        end
    end
    
    % Thermo properties downstream of injector
    s_liq = thermoSat(T_tank(t), 'T', 's_liq');
    s_liq_2 = thermoSat(p_cc(t), 'p', 's_liq');
    s_vap_2 = thermoSat(p_cc(t), 'p', 's_vap');
    h_liq_2 = thermoSat(p_cc(t), 'p', 'h_liq');
    h_vap_2 = thermoSat(p_cc(t), 'p', 'h_vap');
    rho_liq_2 = thermoSat(p_cc(t), 'p', 'rho_liq');
    rho_vap_2 = thermoSat(p_cc(t), 'p', 'rho_vap');
    
    x_2(t) = (s_liq - s_liq_2)/(s_vap_2 - s_liq_2);
    rho_2(t) = 1/(x_2(t)/rho_vap_2 + (1 - x_2(t))/rho_liq_2);
    h_2 = x_2(t)*h_vap_2 + (1 - x_2(t))*h_liq_2;
    
    % Calculate mass flow through the injector
    m_dot_SPI(t) = sqrt(2*rho_liq*(p_tank(t)-p_cc(t)));
    m_dot_HEM(t) = rho_2(t)*sqrt(2*(h_liq-h_2));
    m_dot_ox_in(t) = C_inj * (m_dot_SPI(t) + m_dot_HEM(t))/2;
    
    del_m_tank = -m_dot_ox_in(t)*del_time;
    del_U_tank = (-m_dot_ox_in(t)*h_liq + Q_dot_in)*del_time;
    
    % Step forward in time
    m_tank(t+1) = m_tank(t) + del_m_tank;
    U_tank(t+1) = U_tank(t) + del_U_tank;
    T_tank(t+1) = T_tank(t);
    p_tank(t+1) = p_tank(t);
    x_tank(t+1) = x_tank(t);
    p_cc(t+1) = p_cc(t); % Replace this with combustion chamber model
    t = t+1;
end
time_f = del_time*(length(T_tank)-1);
time = (0:del_time:time_f);
m_dot_ox_in(t) = m_dot_ox_in(t-1); % make mass flow same length as time vector
x_2(t) = x_2(t-1);

%% PLOTS
figure(1)
hold on
plot(time, m_dot_ox_in)
title('Oxidizer Mass Flow vs. Time')
xlabel('Time (s)')
ylabel('Ox Mass Flow (kg/s)')
axis([0 inf 0 inf])
hold off

figure(2)
hold on
plot(time, p_tank./1000000)
title('Oxidizer Tank Pressure vs. Time')
xlabel('Time (s)')
ylabel('Pressure (MPa)')
hold off

figure(3)
hold on
plot(time, T_tank)
title('Oxidizer Tank Temperature vs. Time')
xlabel('Time (s)')
ylabel('Temperature (K)')
hold off

figure(4)
hold on
plot(time, x_tank)
title('Oxidizer Tank Vapour Fraction vs. Time')
xlabel('Time (s)')
ylabel('Vapour Mass Fraction')
hold off

figure(5)
hold on
plot(time, x_2)
title('Vapour Fraction Downstream of Injector (HEM) vs. Time')
xlabel('Time (s)')
ylabel('Vapour Mass Fraction')
hold off

figure(6)
hold on
plot(time, m_tank)
title('Tank Oxidizer Mass vs. Time')
xlabel('Time (s)')
ylabel('Oxidizer Mass (kg)')
hold off

%% SUBFUNCTIONS

    function [varO] = thermoSat(varI,varIname,varOname)
        %THERMOSAT returns thermodynamic properties at saturation

        % Locate which variable column is input, and which is output
        ii = find( (varIname == N2Osat.meta(1,:)), 1, 'first');
        jj = find( (varOname == N2Osat.meta(1,:)), 1, 'first');

        % Interpolate between the two values
        varO = interp1(N2Osat.data(:,ii), N2Osat.data(:,jj), varI);
    end

end



