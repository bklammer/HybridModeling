%% Parametric Study %%
%   See how different parameters change different outputs!

trade = 6; % Which parameter to vary (by index)
n = 20; % How many points to evaluate

% Design parameters are as follows:
var_name = {'Tank Volume (L)', 'Effective Injector Area (mm^2)', ...
    'Fuel Grain Length (cm)', 'Fuel Grain Initial Port Diameter (mm)', ...
    'Nozzle Throat Diameter (mm)', 'Nozzle Area Ratio'};
unit_scale = [1e3, 1e6, 1e2, 1e3, 1e3, 1]; 
lb = [0.004, 1e-5, 0.2, 0.03, 0.02, 3]; % Hardcoded lower bounds
ub = [0.014, 4e-5, 0.6, 0.08, 0.04, 10]; % Hardcoded upper bounds

x = (lb+ub)./2; % Initial point is average of lower and upper bounds
x(7) = 1; % Add seventh element for parametric study identification in Surrocket
x_pert = ones(n,1).*x; % Expand initial point in second dimension

% Vary parameter under examination linearly from lower bnd to upper bnd
x_pert(:,trade) = linspace(lb(trade), ub(trade), n);

len = size(x_pert);
for k = 1:len(1) % For each value
    [~, state] = Surrocket(x_pert(k,:));
    I_sp(k) = state.I_sp;
    max_alt(k) = state.apogee;
    pcc_max(k) = state.pcc_max;
    Thrust_max(k) = state.Thrust_max;
end    

%% PLOT %%

figure
hold on
plot(x_pert(:,trade)*unit_scale(trade), I_sp)
title('Parametric Study - Specific Impulse')
xlabel(var_name{trade})
ylabel('Specific Impulse (s)')
% saveas(gcf, ['Output\Parametric_Study\Isp_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
hold off


figure
hold on
plot(x_pert(:,trade)*unit_scale(trade), max_alt)
title('Parametric Study - Maximum Altitude')
xlabel(var_name{trade})
ylabel('Altitude (m)')
% saveas(gcf, ['Output\Parametric_Study\alt_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
hold off


figure
hold on
plot(x_pert(:,trade)*unit_scale(trade), pcc_max*1e-6)
title('Parametric Study - Maximum Pressure')
xlabel(var_name{trade})
ylabel('Combustion Chamber Pressure (MPa)')
% saveas(gcf, ['Output\Parametric_Study\pcc_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
hold off


figure
hold on
plot(x_pert(:,trade)*unit_scale(trade), Thrust_max)
title('Parametric Study - Maximum Thrust')
xlabel(var_name{trade})
ylabel('Thrust (N)')
% saveas(gcf, ['Output\Parametric_Study\Thrust_AR_', datestr(now, 'yyyy-mm-dd'), '.png']);
hold off