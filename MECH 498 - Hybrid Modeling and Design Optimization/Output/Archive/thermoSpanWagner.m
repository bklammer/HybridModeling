function [out] = thermoSpanWagner(rho, T, in)
%THERMOSPANWAGNER Calculates thermodynamic properties for N2O as a non-ideal gas
%   'in' is a cell array containing the desired output parameters
%   'rho' is the density in kg/m^3
%   'T' is the temperature in K
%   'out' is an array containing the output values in the order listed in 'in'

% Hardcode in data for N2O (from "Modeling Feed System Flow Physics for Self-Pressurizing Propellants)
R = 8.3144598/44.0128*1000; % (J/kg*K)
T_c = 309.52; % (K)
rho_c = 452.0115; % (kg/m^3)
n = [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, ...
    -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038];
n1 = n(1:5); n2 = n(6:12);
a1 = 10.7927224829;
a2 = -8.2418318753;
c0 = 3.5;
v = [2.1769, 1.6145, 0.48393];
u = [879, 2372, 5447];
t = [0.25, 1.125, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5];
d = [1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2];
P = [1, 1, 1, 2, 2, 2, 3];
t1 = t(1:5); t2 = t(6:12);
d1 = d(1:5); d2 = d(6:12);

% Calculate non-dimensional variables
tau = T_c/T;
delta = rho/rho_c;

% Calculate explicit helmholtz energy and derivatives (from 
ao = a1 + a2*tau + log(delta) + (c0-1)*log(tau) + sum(v.*log(1-exp(-u.*tau./T_c)));
ar = sum(n1.*tau.^t1.*delta.^d1) + sum(n2.*tau.^t2.*delta.^d2.*exp(-delta.^P));
ao_tau = a2 + (c0-1)/tau + sum(v.*u./T_c.*exp(-u.*tau./T_c)./(1-exp(-u.*tau./T_c)));
ao_tautau = -(c0-1)/tau.^2 + sum(-v.*u.^2./T_c.^2.*exp(-u.*tau./T_c)./(1-exp(-u.*tau./T_c)).^2);
ar_tau = sum(n1.*t1.*tau.^(t1-1).*delta.^d1) + sum(n2.*t2.*tau.^(t2-1).*delta.^d2.*exp(-delta.^P));
ar_tautau = sum(n1.*t1.*(t1-1).*tau.^(t1-2).*delta.^d1) + sum(n2.*t2.*(t2-2).*tau.^(t2-2).*delta.^d2.*exp(-delta.^P));
ar_delta = sum(n1.*d1.*delta.^(d1-1).*tau.^t1) + sum(n2.*tau.^t2.*delta.^(d2-1).*(d2-P.*delta.^P).*exp(-delta.^P));
ar_deltadelta = sum(n1.*d1.*(d1-1).*delta.^(d1-2).*tau.^t1) + sum(n2.*tau.^t2.*delta.^(d2-2).*((d2-P.*delta.^P).*(d2-1-P.*delta.^P)-P.^2.*delta.^P).*exp(-delta.^P));
ar_deltatau = sum(n1.*d1.*t1.*delta.^(d1-1).*tau.^(t1-1)) + sum(n2.*t2.*tau.^(t2-1).*delta.^(d2-1).*(d2-P.*delta.^P).*exp(-delta.^P));

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

