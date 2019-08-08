function [x] = secant(fun, x1, in)
%SECANT is a zero-finding function based on the secant method
%   'fun' is the function handle for which the zero is desired
%   'x1' is the initial guess
%   'in' is a struct containing any additional inputs 'fun' might require
%   'x' is the value of the zero

x_eps = x1*0.005; % Set the tolerance to be 0.5% of initial guess
x2 = x1-x1*0.01; % Set a second point 1% away from the original guess
F1 = fun(x1, in); % Evaluate function at x1
F2 = fun(x2, in); % Evaluate function at x2
kk = 1; % Set up counter
kk_max = 1000;

while abs(x2-x1)>=x_eps && kk<kk_max % While error is too large and counter is less than max
    x3 = x2 - F2*(x2-x1)/(F2-F1);
    x1 = x2; % Move everything forward
    x2 = x3;
    F1 = F2;
    F2 = fun(x2, in);
    kk = kk+1;
end
x = x2;

end