%% SOBOL' SENSITIVITY ANALYSIS %%
% Coded according to the algorithm provided in chapter 4 of "Global
% Sensitivity Analysis: A Primer" for variance based methods

tic;

N = 1280; % Number of runs
k = 4; % Dimensionality of uncertainty vector
n_eval = N*(k+2); % Total number of model evaulations

S0 = 1000;
fun = @(x) x(1)*x(2)*S0-x(3)-x(4); % Specify which function to quantify sensitivity
%ub = ones(N,k).*[1000, 100, 10]; % Upper bounds
%lb = ones(N,k).*[1, 1, 1]; % lower bounds

s = sobolset(2*k); % Define a Sobol' set
A = s(2:(N+1),1:k); % Split the set into two equally-sized matrices
B = s(2:(N+1),(k+1):(2*k)); % Split the set into two equally-sized matrices

A(:,2) = icdf('beta',A(:,2),0.2,15); % Replace with inverse cumulative beta distribution (assume sobol points are outputs of cdf)
B(:,2) = icdf('beta',B(:,2),0.2,15);

%A = A.*(ub-lb)+lb; % Apply upper and lower bounds, assuming uniform distribution
%B = B.*(ub-lb)+lb; % Apply upper and lower bounds, assuming uniform distribution

C = zeros(N,k,k);
for i=1:k
    C(:,:,i) = B;
    C(:,i,i) = A(:,i);
end

yA = zeros(N,1);
yB = zeros(N,1);
yC = zeros(N,k);
for j=1:N
    yA(j) = fun(A(j,:));
    yB(j) = fun(B(j,:));
    for i=1:k
        yC(j,i) = fun(C(j,:,i));
    end
end

f0 = mean(yA);
for i=1:k
    S(i) = (mean(yA.*yC(:,i))-f0^2)./(mean(yA.^2)-f0^2); % First order sensitivity indices
    ST(i) = 1-(mean(yB.*yC(:,i))-f0^2)./(mean(yA.^2)-f0^2); % Total effect indices
end
    
toc;