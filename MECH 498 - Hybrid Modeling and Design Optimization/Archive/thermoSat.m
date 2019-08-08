function [varO] = thermoSat(varI,varIname,varOname)
%THERMOSAT returns thermodynamic properties at saturation

N2Osat = load('N2Osat.mat');
N2Osat = N2Osat.N2Osat;

% Locate which variable column is input, and which is output
ii = find( (varIname == N2Osat.meta(1,:)), 1, 'first');
jj = find( (varOname == N2Osat.meta(1,:)), 1, 'first');

% Interpolate between the two values
varO = interp1(N2Osat.data(:,ii), N2Osat.data(:,jj), varI);

end

