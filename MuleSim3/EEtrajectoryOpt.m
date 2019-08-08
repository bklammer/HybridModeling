function [D] = EEtrajectoryOpt(x)
%EETRAJECTORYOPT Calculates the distance between a set of generated trajectories

% Load trajectories?

% Find distance between trajectories
d = zeros(M);
for m=1:M
    for l=(m+1):M
        for i=1:(k+1) % Runs through points in mth trajectory
            for j=1:(k+1) % Runs through points in lth trajectory
                temp(i,j) = sqrt(sum((X(i,:,m)-X(j,:,l)).^2)); % Sum of distance between points
            end
        end
        d(m,l) = sum(temp, 'all'); % Distance between trajectories m and l
    end
end
d = d+d'; % Add transpose to make indexing more reliable;


% Find combination that maximizes distance between trajectories


n = nchoosek(1:M,r); % All possible combinations
D_trajectory = zeros(length(n),1);
for i=1:length(n)
    ind = nchoosek(n(i,:),2); % Find indices
    dtemp = diag(d(ind(:,1),ind(:,2))); % Find distances between individual specified trajectoris
    D_trajectory(i) = sqrt(sum(dtemp.^2)); % Calculate total distance between trajectories
end
[~,ind] = max(D_trajectory);
T = zeros((k+1),k,r);
for i=1:r
    T(:,:,i) = X(:,:,n(ind,i));
end


end

