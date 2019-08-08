function [p] = CombustionChamber(p)
%COMBUSTIONCHAMBER Models a hybrid rocket combustion chamber
%   'p' is a struct containing at least the following fields:
%       t, r_cc, m_dot_ox_in, a, n, zeta_d, A_th, T_stag, R_cc, k_cc, T_cc,
%       zeta_cstar, del_time, r_dot, m_dot_f, m_dot_cc, OF, p_stag, p_cc,
%       vel_cc, m_ox_tank, t_max, rho_f

% Assign relevant fields in struct to local variables for code legibility, 
% accounting, and debugging purposes
r_cc = p.r_cc(p.t);
m_dot_ox_in = p.m_dot_ox_in(p.t);
a = p.a;
n = p.n;
rho_f = p.rho_f;
L = p.L;
zeta_d = p.zeta_d; % (-) Discharge correction factor
A_th = p.A_th;
T_stag = p.T_stag(p.t);
R_cc = p.R_cc(p.t);
k_cc = p.k_cc(p.t);
T_cc = p.T_cc(p.t);
zeta_cstar = p.zeta_cstar; % (-) Characteristic velocity correction factor
del_time = p.del_time; % (s) Simulation time step

%% Combustion Chamber Calculations
A_cc = pi()*r_cc^2; % Port geometry
del_r_cc = a*(m_dot_ox_in/A_cc)^n * del_time; % Regression rate law
m_dot_f = 2*pi()*r_cc*L*rho_f*(del_r_cc/del_time); % Fuel mass flow rate from geometry

m_dot_cc = m_dot_f + m_dot_ox_in; % Total mass flow in
OF = m_dot_ox_in/m_dot_f; % Calculate Oxidizer-Fuel Ratio
p_stag = m_dot_cc/(zeta_d*A_th)*sqrt(T_stag*R_cc/k_cc*((k_cc+1)/2)^((k_cc+1)/(k_cc-1))); % Steady state choked flow expression for pressure
p_cc = p_stag*((T_cc/T_stag).^(k_cc/(k_cc-1)));

thermo_comb = CEAProp(OF, p_cc, {'T', 'rho', 'cp', 'k', 'M'});
thermo_comb = num2cell(thermo_comb);
[T_stag, rho_cc, cp_cc, k_cc, M_cc] = thermo_comb{:};
R_cc = 8.3144598/M_cc;

vel_cc = m_dot_cc/(rho_cc*A_cc);

T_stag = T_stag*zeta_cstar^2;
T_cc = T_stag - vel_cc^2/(2*cp_cc);
p_stag = p_cc*(T_stag/T_cc).^(k_cc/(k_cc-1)); % Re-calculate p_stag using new T_cc, T_stag, and k_cc


%% Assignment to struct and time step
% Assign values of interest to output struct
p.r_dot(p.t) = del_r_cc/del_time;
p.m_dot_f(p.t) = m_dot_f;
p.m_dot_cc(p.t) = m_dot_cc;
p.OF(p.t) = OF;
p.p_stag(p.t) = p_stag;
p.p_cc(p.t) = p_cc;
p.R_cc(p.t) = R_cc;
p.T_stag(p.t) = T_stag;
p.vel_cc(p.t) = vel_cc;
p.T_cc(p.t) = T_cc;

if p.m_ox_tank(p.t)/p.m_ox_tank(1) > 0.03 && p.t < p.t_max
% If less than max burn time and more than 3% of oxidizer is left, step forward in time
    % Assign values at next time step to ouptut struct
    p.r_cc(p.t+1) = r_cc + del_r_cc;
    p.p_cc(p.t+1) = p_cc;
    p.T_cc(p.t+1) = T_cc;
    p.T_stag(p.t+1) = T_stag;
    p.R_cc(p.t+1) = R_cc;
    p.k_cc(p.t+1) = k_cc;
end


%% Subfunctions
function [out] = CEAProp(OF, pressure, in)
%MASSFRAC Calculates the mass fractions of products for a given reaction
%   'OF' is the oxidizer to fuel weight ratio
%   'p' is the pressure of the combustion chamber
%   'alpha' is an array of weight fractions
    if OF < p.CEA.OF(1) % Check that query is inside bounds
        error('CEAProp:lowOF', 'Outside O/F range: \n OF = %f', OF)
    elseif OF > p.CEA.OF(end)
        error('CEAProp:highOF', 'Outside O/F range: \n OF = %f', OF)
    elseif pressure < p.CEA.p(1)
        error('CEAProp:lowP', 'Outside pressure range: \n p = %f Pa', pressure)
    elseif pressure > p.CEA.p(end)
        error('CEAProp:highP', 'Outside pressure range: \n p = %f Pa', pressure)
    end
    for mm = 1:length(p.CEA.OF)
        if OF < p.CEA.OF(mm) % Find first index that is larger than OF
            break
        end
    end
    for nn = 1:length(p.CEA.p)
        if pressure < p.CEA.p(nn) % Find first index that is larger than p
            break
        end
    end
    out = zeros(size(in));
    for kk = 1:length(in) % Iterate through all products
        % Do some bilinear interpolation (equations from wikipedia)
        OF_int = [p.CEA.OF(mm)-OF, OF-p.CEA.OF(mm-1)];
        p_int = [p.CEA.p(nn)-pressure; pressure-p.CEA.p(nn-1)];
        C_int = 1/((p.CEA.OF(mm)-p.CEA.OF(mm-1))*(p.CEA.p(nn)-p.CEA.p(nn-1)));
        out_int = p.CEA.(in{kk})(mm-1:mm,nn-1:nn);
        out(kk) = C_int.* OF_int*out_int*p_int;
    end
end % end CEAProp

end