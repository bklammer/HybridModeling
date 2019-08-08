function [MAT] = CEAtoMATLAB(folder, file, thermo)
%CEAtoMATLAB takes a CEA output file and converts it to a .mat file
% 'folder' is a string specifying a folder name to look for file in
% 'file' is a desired file name for the saved .mat file
% 'thermo' is a CEA thermodynamic database saved as a .mat file

% Handle different inputs
if nargin == 0
    folder = '';
    file = 'MuleSim2CEA';
    thermo = load('MuleSimTHERMO.mat'); thermo = thermo.MAT;
elseif nargin == 1
    file = 'MuleSim2CEA';
    thermo = load('MuleSimTHERMO.mat'); thermo = thermo.MAT;
elseif nargin == 2
    thermo = load('MuleSimTHERMO.mat'); thermo = thermo.MAT;
end

% Import CEA '.out' file to string cell array
[fileName, filePath] = uigetfile([folder, '*.out']);
if isequal(fileName, 0)
    error('User selected cancel')
end
fullFileName = [filePath, fileName];
fileID = fopen(fullFileName);
C = textscan(fileID, '%s');
CEA = C{1};

% Loop initializations
MAT = struct();
MAT.OF(1) = 0;
k = 1;
m = 1;
n = 1;
O = 1;
F = 1;
while k < length(CEA)

    if O == 1 && strcmp(CEA{k}, 'OXIDANT') % Extract oxidizer data
        k=k+1;
        while ~strcmp(CEA{k}, 'FUEL') && k < length(CEA)
            oxName = CEA{k};
            oxName = strrep(oxName, '+', '_plus');
            oxName = strrep(oxName, '-', '_minus');
            oxName = matlab.lang.makeValidName(oxName);
            k=k+1;
            MAT.OX.(oxName).wtfrac = str2double(CEA{k});
            k=k+1;
            MAT.OX.(oxName).h = str2double(CEA{k});
            k=k+1;
            MAT.OX.(oxName).T = str2double(CEA{k});
            k=k+1;
        end
        O = 0;
    end
    
    if F == 1 && strcmp(CEA{k}, 'FUEL') % Extract fuel data
        k=k+1;
        while ~strcmp(CEA{k}, 'O/F=') && k < length(CEA)
            fuelName = CEA{k};
            fuelName = strrep(fuelName, '+', '_plus');
            fuelName = strrep(fuelName, '-', '_minus');
            fuelName = matlab.lang.makeValidName(fuelName);
            k=k+1;
            MAT.FUEL.(fuelName).wtfrac = str2double(CEA{k});
            k=k+1;
            MAT.FUEL.(fuelName).h = str2double(CEA{k});
            k=k+1;
            MAT.FUEL.(fuelName).T = str2double(CEA{k});
            k=k+1;
        end
        F = 0;
    end
    
    if strcmp(CEA{k}, 'O/F=') % Store O/F data
        k=k+1;
        OF = str2double(CEA{k});
        if MAT.OF(m) ~= OF && MAT.OF(m) ~= 0
            m=m+1;
            MAT.OF(m) = OF;
            n = 1;
        else
            MAT.OF(m) = OF;
        end
    end
    
    if m == 1 && strcmp(CEA{k}, 'BAR') % Store pressure data
        k=k+1;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.P(n_temp) = str2double(CEA{k});
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'FRACTIONS') % Extract product mass fractions
        k=k+1;
        while k < length(CEA) && ~strcmp(CEA{k}, '*')
            if isnan(str2double(CEA{k}))
                prodName = CEA{k};
                prodName = strrep(prodName, '+', '_plus');
                prodName = strrep(prodName, '-', '_minus');
                prodName = erase(prodName, '*');
                prodName = matlab.lang.makeValidName(prodName);
                n_temp = n;
                k=k+1;
            else
                MAT.PRODUCTS.(prodName).MASSFRAC(m,n_temp) = str2double(CEA{k});
                k=k+1;
                n_temp=n_temp+1;
            end
        end
        n = n_temp;
    end
    
    k=k+1;
end

% Ensure all reactant mass fraction matrices are the same size
m = length(MAT.OF);
n = length(MAT.P);
prodNames = fieldnames(MAT.PRODUCTS);
for k = 1:length(prodNames)
    [mp,np] = size(MAT.PRODUCTS.(prodNames{k}).MASSFRAC);
    if mp < m || np < n
        MAT.PRODUCTS.(prodNames{k}).MASSFRAC(m,n) = 0;
    end
end

% Generate thermodynamic property lookup table for products
for k = 1:length(prodNames) % Iterate through each product species
    tp = thermo.(prodNames{k}); % Set up temporary variable to improve legibility
    M = tp.prop(3)/1000; % (kg/mol)
    Rk = 8.3144598/M; % (J/kg*K)
    
    T = tp.T(1):50:tp.T(end); % Change the middle number to define the fineness of the temperature step
    
    ii = 1;
    mm = zeros(size(T));
    for kk = 1:length(T) % Sets mm (index for temperature range) to correct range for range T
        if T(kk) <= tp.T(ii,2)
            mm(kk) = ii;
        else
            ii = ii + 1;
            mm(kk) = ii;
        end
    end
    
    cp = zeros(size(T));
    h = zeros(size(T));
    for kk = 1:length(T)
        cp(kk) = Rk*sum( tp.a(mm(kk),:) .* (T(kk) .^ (tp.n(mm(kk),:))) ); % (J/kg*K)
        h_temp = (tp.a(mm(kk),:)./(tp.n(mm(kk),:)+1)) .* (T(kk) .^ (tp.n(mm(kk),:)+1)); % Uses the formula for integration of a polynomial, since delh = integral(cp)dT
        h_temp(~isfinite(h_temp)) = tp.a(mm(kk), ~isfinite(h_temp)).* log(T(kk)); % Replaces the infinite term with ln(T) since polynomial integration formula doesn't work for 1/T
        h(kk) = Rk*(sum(h_temp) + tp.b(mm(kk),1)); % (J/kg)
    end
    
    MAT.PRODUCTS.(prodNames{k}).T = T;
    MAT.PRODUCTS.(prodNames{k}).cp = cp;
    MAT.PRODUCTS.(prodNames{k}).h = h;
    MAT.PRODUCTS.(prodNames{k}).M = tp.prop(3);
end

% Add molecular weight to oxidizer
oxNames = fieldnames(MAT.OX);
for k = 1:length(oxNames)
    MAT.OX.(oxNames{k}).M = thermo.(oxNames{k}).prop(3);
end

% Add molecular weight to fuel
fuelNames = fieldnames(MAT.FUEL);
for k = 1:length(oxNames)
    MAT.FUEL.(fuelNames{k}).M = thermo.(fuelNames{k}).prop(3);
end

% Save as '.mat' file
save([folder file '.mat'],'MAT')
fclose(fileID);
end

