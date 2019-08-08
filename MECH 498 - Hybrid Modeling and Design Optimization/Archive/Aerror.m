function A = Aerror(M, in) % Finds the difference between the estimated and actual nozzle ratio
    A = (1./M.^2).*(2./(in.k+1).*(1+(in.k-1)./2.*M.^2)).^((in.k+1)./(in.k-1)) - in.A^2;
end