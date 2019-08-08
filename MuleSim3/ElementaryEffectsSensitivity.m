%% ELEMENTARY EFFECTS SENSITIVITY ANALYSIS %%
% Coded according to the algorithms provided in chapter 3 of "Global
% Sensitivity Analysis: A Primer" for elementary effects methods

clear all;
load('ElementaryEffectsTrajectory_2019-04-10.mat');

tic;

% Order is as follows:
%   V_tank, tank_level, p_tank_init, C_inj, L, d_port_init, d_th,
%   A_ratio_nozzle, p_feed, rho_f, a, n, launchAlt, p_amb, T_amb, zeta_d,
%   zeta_cstar, zeta_CF, C_D
average = [0.009366004, 0.6, 5000000, 0.0000213, 0.251541, 0.062527, 0.025552, ...
    5.986583, 100000, 930, 0.000155, 0.5, 1400, 85600, 301, 1.05, 0.818, 0.9, 1];
deviation = [0.02, 0.10, 0.05, 0.05, 0.01, 0.03, 0.07, 0.14, 0.5, 0.03, 0.10, ...
    0.01, 0.007, 0.10, 0.02, 0.05, 0.10, 0.05, 0.10];


% % Specify which function to quantify sensitivity
% a = [78, 12, 0.5, 2, 97, 33];
% fun = @(x) prod( (abs(4.*x-2)+a)./(1+a) );
fun = @Surrocket;

% Trajectory definition
k = 19; % Dimensionality of uncertainty vector
delta = 2/3; % Change per variable
r = 8; % Desired number of trajectories
M = 25; % Number of generated trajectories

% % Build trajectories
% J = ones((k+1),k); % Generate matrix of ones
% B = tril(J,-1); % Generate lower triangular matrix of ones
% for m=1:M
%     x = (unidrnd(ones(1,k)*2)-1)/3; % Generate random base value vector
%     D = diag((rand(k,1) > 0.5)*2 - 1); % Diagonal matrix of ones or negative ones
%     ind = randperm(k);
%     P = zeros(k);
%     for i=1:k
%         P(i,ind(i)) = 1; % Random permutation matrix
%     end
%     X(:,:,m) = (J.*x + delta/2*((2*B-J)*D + J))*P;
% end
% 
% 
% % Find distance between trajectories
% d = zeros(M);
% for m=1:M
%     for l=(m+1):M
%         for i=1:(k+1) % Runs through points in mth trajectory
%             for j=1:(k+1) % Runs through points in lth trajectory
%                 temp(i,j) = sqrt(sum((X(i,:,m)-X(j,:,l)).^2)); % Sum of distance between points
%             end
%         end
%         d(m,l) = sum(temp, 'all'); % Distance between trajectories m and l
%     end
% end
% d = d+d'; % Add transpose to make indexing more reliable;
% 
% 
% % Find combination that maximizes distance between trajectories
% 
% 
% n = nchoosek(1:M,r); % All possible combinations
% D_trajectory = zeros(length(n),1);
% for i=1:length(n)
%     ind = nchoosek(n(i,:),2); % Find indices
%     dtemp = diag(d(ind(:,1),ind(:,2))); % Find distances between individual specified trajectoris
%     D_trajectory(i) = sqrt(sum(dtemp.^2)); % Calculate total distance between trajectories
% end
% [~,ind] = max(D_trajectory);
% T = zeros((k+1),k,r);
% for i=1:r
%     T(:,:,i) = X(:,:,n(ind,i));
% end




% Convert uniform distribution to appropriate distribution
V = T;
V(:,1,:) = T(:,1,:)*0.6826 + 0.1587; % V_tank, normal distribution
V(:,2,:) = T(:,2,:)*0.6826 + 0.1587; % tank_level, normal distribution
V(:,5,:) = T(:,5,:)*0.6826 + 0.1587; % L, normal distribution
V(:,6,:) = T(:,6,:)*0.6826 + 0.1587; % d_port_init, normal distribution
V(:,7,:) = T(:,7,:)*0.6826 + 0.1587; % d_th, normal distribution
V(:,8,:) = T(:,8,:)*0.6826 + 0.1587; % A_ratio_nozzle, normal distribution
V(:,13,:) = T(:,13,:)*0.6826 + 0.1587; % launchAlt, normal distribution
V(:,14,:) = T(:,14,:)*0.6826 + 0.1587; % p_amb, normal distribution
V(:,15,:) = T(:,15,:)*0.6826 + 0.1587; % T_amb, normal distribution
% All others are uniform distributions

V(:,1,:) = icdf('normal',V(:,1,:),0,1); % V_tank, normal distribution
V(:,2,:) = icdf('normal',V(:,2,:),0,1); % tank_level, normal distribution
V(:,3,:) = 2*V(:,3,:)-1; % p_tank_init, uniform distribution
V(:,4,:) = 2*V(:,4,:)-1; % C_inj, uniform distribution
V(:,5,:) = icdf('normal',V(:,5,:),0,1); % L, normal distribution
V(:,6,:) = icdf('normal',V(:,6,:),0,1); % d_port_init, normal distribution
V(:,7,:) = icdf('normal',V(:,7,:),0,1); % d_th, normal distribution
V(:,8,:) = icdf('normal',V(:,8,:),0,1); % A_ratio_nozzle, normal distribution
V(:,9,:) = 2*V(:,9,:)-1; % p_feed, uniform distribution
V(:,10,:) = 2*V(:,10,:)-1; % rho_f, uniform distribution
V(:,11,:) = 2*V(:,11,:)-1; % a, uniform distribution
V(:,12,:) = 2*V(:,12,:)-1; % n, uniform distribution
V(:,13,:) = icdf('normal',V(:,13,:),0,1); % launchAlt, normal distribution
V(:,14,:) = icdf('normal',V(:,14,:),0,1); % p_amb, normal distribution
V(:,15,:) = icdf('normal',V(:,15,:),0,1); % T_amb, normal distribution
V(:,16,:) = 2*V(:,16,:)-1; % zeta_d, uniform distribution
V(:,17,:) = 2*V(:,17,:)-1; % zeta_cstar, uniform distribution
V(:,18,:) = 2*V(:,18,:)-1; % zeta_CF, uniform distribution
V(:,19,:) = 2*V(:,19,:)-1; % C_D, uniform distribution

% Convert non-dimensional values to dimensional values
V = average.*(deviation.*V + 1);

% Evaluate model at specified points
for i=1:r
    for j=1:(k+1)
        Y(j,i) = fun(V(j,:,i));
    end
end

% Calculate elementary effects
EE = zeros(r,k);
for i=1:r
    for j=2:(k+1)
        EEtemp = (Y(j,i)-Y(j-1,i))./(T(j,:,i)-T(j-1,:,i));
        ind = ~isinf(EEtemp);
        EE(i,ind) = EEtemp(ind);
    end
end

% Calculate sensitivity measures
mu = (nanmean(EE))';
sigma = (nanstd(EE))';
mu_star = (nanmean(abs(EE)))';

results = [mu_star, mu, sigma];

% toc;