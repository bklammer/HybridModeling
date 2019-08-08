fval_sz = ((median(fval_vel)-fval_vel)*200+50)';
fval_sz(fval_sz<30) = 30;

x_scaled = x_vel./ub;

one = ones(size(fval_sz));
figure('Position', [100 100 400 300])
hold on;
scatter([1,2,3,4,5,6],x_ga_scaled,50,'o')

for i=1:6
    scatter((i.*one), x_scaled(:,i),fval_sz,'.')
end

xlim([0,7]);
ylim([0,1]);
title('Pattern Search Results');
designVars = {'Vtank', 'Cinj', 'L', 'dport', 'dth', 'AR'};
set(gca, 'XTick', 1:6, 'XTickLabel', designVars);