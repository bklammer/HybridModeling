function [rho_ox_2, T_ox_2] = multiVariateNewtonRaphson(rho_tank, rho_ox_2, T_tank, T_ox_2, s, p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x2(1) = rho_tank;
x2(2) = T_tank;
x1(1) = rho_ox_2;
x1(2) = T_ox_2;

F1 = inj_2(x1);
F2 = inj_2(x2);
J = zeros(2);
for ii=1:2 % Construct the Jacobian using secant method
    for jj=1:2
        x_pert = x2;
        x_pert(jj) = x1(jj);
        F_pert = inj_2(x_pert);
        J(ii,jj) = (F2(ii)-F_pert(ii)) / (x2(jj)-x1(jj));
    end
end
Jinv = J^-1;
n = 1;
n_max = 1000;
x_eps = [1, 1];
while all(abs(x2-x1)>x_eps) && n<n_max
    x3 = x2 - (Jinv*F2')'; % Calculate new iteration
    x1 = x2; % Move everything forward
    x2 = x3;
    F1 = F2;
    F2 = inj_2(x2);
    Jinv = Jinv + ((x2-x1)'-Jinv*(F2-F1)')/((x2-x1)*Jinv*(F2-F1)')*(x2-x1)*Jinv; % Sherman-Morrison formula to update inverse of Jacobian
    n = n+1;
end

rho_ox_2 = x2(1);
T_ox_2 = x2(2);

function F = inj_2(x)
    F(1) = thermoSpanWagner({'s'}, x(1), x(2)) - s;
    F(2) = thermoSpanWagner({'p'}, x(1), x(2)) - p;
end

end

    

