% Converts unevenly spaced occasionally repeating time curves to evenly
% spaced time curves sampeld at 100Hz

n = 2;
while n < length(sample)+1
    if sample(n-1, 1) == sample(n, 1)
        sample(n,:) = [];
        n = n -1;
    end
    n = n+1;
end
t = 0:0.01:sample(end,1);
data = (interp1(sample(:,1), sample(:,2), t))';