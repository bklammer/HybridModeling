function [p] = Nozzle(p)
%NOZZLE Models an ideal rocket nozzle
%   'p' is a struct containing at least the following fields:
%       t, A_ratio_nozzle, Ma, p_stag, k_cc, x_tank, R_cc, T_stag, zeta_CF,
%       m_dot_cc, p_atm, A_th, p_exit, vel_exit, F_thrust

% Assign relevant fields in struct to local variables for code legibility, 
% accounting, and debugging purposes
A_ratio_nozzle = p.A_ratio_nozzle;
A_ratio_nozzle_eps = p.A_ratio_nozzle_eps;
Ma = p.Ma;
p_stag = p.p_stag(p.t);
k_cc = p.k_cc(p.t);
x_tank = p.x_tank(p.t);
R_cc = p.R_cc(p.t);
T_stag = p.T_stag(p.t);
zeta_CF = p.zeta_CF;
m_dot_cc = p.m_dot_cc(p.t);
p_atm = p.p_atm;
A_th = p.A_th;


%% Nozzle and Thrust Calculations
Ain.k = k_cc;
Ain.A = A_ratio_nozzle;
if abs(Aerror(Ma, Ain)) > A_ratio_nozzle_eps % Nozzle mach number solver
    Ma = secant(@Aerror, Ma, Ain);
end
p_exit = p_stag./(1+(k_cc-1)./2.*Ma.^2).^(k_cc/(k_cc-1)); % (Pa) Pressure of gas at exit of nozzle
if x_tank >= 1 % Physically very wrong, but the condition is used to hide the transition that results from overexpansion
    p_exit = p_stag*(2/(k_cc+1))^(k_cc/(k_cc-1)); % Throat pressure
    vel_exit = sqrt(2*k_cc*R_cc*T_stag/(k_cc+1)); % Throat velocity
    A_ratio_nozzle_eff = 1; % Throat area ratio
else
    T_exit = T_stag./(1+(k_cc-1)./2.*Ma.^2); % (K) Temperature of gas at exit of nozzle
    vel_exit = Ma.*sqrt(k_cc.*R_cc.*T_exit); % (m/s) Velocity of gas at exit of nozzle
    A_ratio_nozzle_eff = A_ratio_nozzle; % Effective nozzle area ratio is actual nozzle area ratio
end
F_thrust = zeta_CF*(m_dot_cc.*vel_exit + (p_exit-p_atm).*A_th.*A_ratio_nozzle_eff); % (N) Rocket Motor Thrust!!!


%% Assign values if interest to output struct
p.Ma = Ma;
p.p_exit(p.t) = p_exit;
p.vel_exit(p.t) = vel_exit;
p.F_thrust(p.t) = F_thrust;

end % end of main program


%% Finds the difference between the estimated and actual nozzle ratio
function A = Aerror(M, in) 
    A = (1/M^2)*(2./(in.k+1).*(1+(in.k-1)./2.*M.^2)).^((in.k+1)./(in.k-1)) - in.A^2;
end