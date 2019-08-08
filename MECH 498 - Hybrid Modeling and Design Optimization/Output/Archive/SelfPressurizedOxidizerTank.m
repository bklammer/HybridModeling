function [p] = SelfPressurizedOxidizerTank(p)
%SELFPRESSURIZEDOXIDIZERTANK Models a self-pressurized oxidizer tank
%   'p' is a struct containing at least the following fields:
%       t, x_tank, U_tank, m_ox_tank, V_tank, V_tank_eps, T_tank, p_cc,
%       C_inj, del_time, p_tank, m_dot_ox_in, t_max, N2Osat

% Assign relevant fields in struct to local variables for code legibility, 
% accounting, and debugging purposes
x_tank = p.x_tank(p.t); % (-) Tank vapour mass fraction
U_tank = p.U_tank(p.t); % (J) Tank internal energy 
m_ox_tank = p.m_ox_tank(p.t); % (kg) Tank mass
V_tank = p.V_tank; % (m^3) Tank volume
V_tank_eps = p.V_tank_eps; % (m^3) Tank volume acceptable error
u_tank_eps = p.u_tank_eps;
T_tank = p.T_tank(p.t); % (K) Tank temperature guess
p_cc = p.p_cc(p.t); % (Pa) Combuston chamber pressure
C_inj = p.C_inj; % (m^2) Injector area coefficient
del_time = p.del_time; % (s) Simulation time step

if x_tank < 1
    % OXIDIZER TANK CALCULATIONS FOR SATURATED LIQUID AND VAPOUR %
    Vin.U = U_tank;
    Vin.m = m_ox_tank;
    Vin.V = V_tank;
    if abs(Verror(T_tank, Vin)) > V_tank_eps
        T_tank = secant(@Verror, T_tank, Vin);
    end

    % Calculation of tank thermo properties
    thermo_sat_tank = thermoSat(T_tank, 'T', {'p', 'h_liq', 'rho_liq', 'u_liq', 'u_vap'});
    thermo_sat_tank = num2cell(thermo_sat_tank);  
    [p_tank, h_liq, rho_liq, u_liq, u_vap] = thermo_sat_tank{:};
    x_tank = (U_tank/m_ox_tank - u_liq)/(u_vap - u_liq);

    % Thermo properties downstream of injector for HEM model
    s_liq = thermoSat(T_tank, 'T', {'s_liq'});
    thermo_sat_inj = thermoSat(p_cc, 'p', {'s_liq', 's_vap', 'h_liq', 'h_vap', 'rho_liq', 'rho_vap'});
    thermo_sat_inj = num2cell(thermo_sat_inj);
    [s_liq_2, s_vap_2, h_liq_2, h_vap_2, rho_liq_2, rho_vap_2] = thermo_sat_inj{:};
    x_2 = (s_liq - s_liq_2)/(s_vap_2 - s_liq_2); % Mass fraction downstream assuming isentropic expansion at saturation (only use for HEM mass flow)
    rho_ox_2 = 1/(x_2/rho_vap_2 + (1 - x_2)/rho_liq_2);
    h_2 = x_2*h_vap_2 + (1 - x_2)*h_liq_2;

    % Calculate mass flow through the injector
    G_SPI = sqrt(2*rho_liq*(p_tank-p_cc));
    G_HEM = rho_ox_2*sqrt(2*(h_liq-h_2));
    m_dot_ox_in = C_inj*(G_SPI+G_HEM)/2;

    del_m_ox_tank = -m_dot_ox_in*del_time;
    del_U_tank = -m_dot_ox_in*h_liq*del_time;

    if x_tank > 1
        p.burn_time = p.t*p.del_time; % Record burn time to be the time at which the tank runs out of oxidizer
    end
else
    % OXIDIZER TANK CALCULATIONS FOR VAPOUR ONLY %
    rho_tank = m_ox_tank./V_tank;
    u_tank = U_tank/m_ox_tank;
    % Assign variables to struct for input to 'uerror'
    uin.rho = rho_tank;
    uin.u = u_tank;
    if abs(uerror(T_tank, uin)) > u_tank_eps
        T_tank = secant(@uerror, T_tank, uin);
    end

    thermo_span_tank = thermoSpanWagner(rho_tank, T_tank, {'p','h'});
    p_tank = thermo_span_tank(1);
    h_tank = thermo_span_tank(2) + 7.3397e+05; % Convert from Span-Wagner enthalpy convention to NIST
    m_dot_ox_in = C_inj*sqrt(2*rho_tank*(p_tank-p_cc)); % Incompressible fluid assumption (better than nothing)

    del_m_ox_tank = -m_dot_ox_in*del_time;
    del_U_tank = -m_dot_ox_in*h_tank*del_time;
end


% Assign values if interest to output struct
p.p_tank(p.t) = p_tank;
p.m_dot_ox_in(p.t) = m_dot_ox_in;

if p.m_ox_tank(p.t)/p.m_ox_tank(1) > 0.03 && p.t < p.t_max
% If less than max burn time and more than 3% of oxidizer is left, step forward in time
    % Assign values at next time step to ouptut struct
    p.m_ox_tank(p.t+1) = m_ox_tank + del_m_ox_tank;
    p.U_tank(p.t+1) = U_tank + del_U_tank;
    p.x_tank(p.t+1) = x_tank;
    p.T_tank(p.t+1) = T_tank;
end


%% Subfunctions
function [val_out] = thermoSat(val_in, in, out)
%THERMOSAT returns thermodynamic properties at saturation
%   'val_in' is a double specifying the value of thermodynamic property inputted
%   'in' is a string specifying the given input thermodynamic property 
%   'out' as a cell array specifying the desired output thermodynamic property
%   'val_out' is a double specifying the value of thermodynamic property returned
    val_out = zeros(size(out));
    for mm = 1:length(out) % For each desired output
        % Locate which variable column is input, and which is output
        ii = find( (in == p.N2Osat.meta(1,:)), 1, 'first');
        jj = find( (out(mm) == p.N2Osat.meta(1,:)), 1, 'first');
        % Interpolate between the two values
        for kk = 1:length(p.N2Osat.data(:,ii))
            if val_in < p.N2Osat.data(kk,ii)
                break
            end
        end
        val_out(mm) = (val_in - p.N2Osat.data(kk-1,ii))./(p.N2Osat.data(kk,ii)-p.N2Osat.data(kk-1,ii)).*(p.N2Osat.data(kk,jj)-p.N2Osat.data(kk-1,jj)) + p.N2Osat.data(kk-1,jj);
    end
end % end thermoSat

%% Thermodynamic Properties of N2O as a Real Gas
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


%% Finds the difference between the estimated and actual tank volume
function V = Verror(T, in) 
    thermo = thermoSat(T, 'T', {'rho_liq','rho_vap','u_liq','u_vap'});
    thermo = num2cell(thermo);
    [rho_l,rho_v,u_l,u_v] = thermo{:};
    x = (in.U/in.m - u_l)/(u_v - u_l);
    V = in.m*((1-x)/rho_l + x/rho_v) - in.V;
end


%% Finds the difference between the estimated and actual tank internal energy
function U = uerror(T, in) 
    U = (thermoSpanWagner(in.rho, T, {'u'}) + 7.3397e+5) - in.u;
end

end % end loop






