function [MAT] = CEAtoMATLAB2(file, folder)
%CEAtoMATLAB takes a CEA output file and converts it to a .mat file
% 'folder' is a string specifying a folder name to look for file in
% 'file' is a desired file name for the saved .mat file

% Handle different inputs
if nargin == 0
    folder = '';
    file = 'MuleSim2CEA';
elseif nargin == 1
    folder = '';
end

%% READING COMBUSTION PROPERTIES FILE
% Import CEA '.out' file to string cell array
[fileName, filePath] = uigetfile([folder, '*.out'], 'Please select a CEA output file');
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
            if strcmp(CEA{k}, 'OXIDANT')
                k=k+1;
            end
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
            if strcmp(CEA{k}, 'FUEL')
                k=k+1;
            end
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

%% READING THERMODYNAMIC DATABASE FILE
% Import 'thermo.inp' file to char array
[fileName, filePath] = uigetfile([folder, '*.inp'], 'Please select a CEA thermodynamic database');
if isequal(fileName, 0)
    error('User selected cancel')
end
fullFileName = [filePath, fileName];
C = regexp(fileread(fullFileName), '\r?\n', 'split');
CEA = char(C{:});
CEA = CEA(:,1:80);

% Loop through CEA looking for products
k = 1;
while k < length(CEA)
    b1 = ~strcmp(CEA(k,1), {' ', '!', '#', '-'});
    b2 = ~strcmp(CEA(k,1:18), {'thermo            ', 'END PRODUCTS      ', 'END REACTANTS     '});
    
    if all(b1) && all(b2) % If the line isn't a comment, a space, or a keyword
        name = CEA(k,1:18); % Extract name of species
        name = strrep(name, '+', '_plus');
        name = strrep(name, '-', '_minus');
        name = matlab.lang.makeValidName(name);
%         if length(name) > 15
%             name = name(1:15); % Limit it to 15 characters because the CEA output is limited to 15 characters
%         end
        k=k+1;
        
        if any(strcmp(name, prodNames))
            numT = str2double(CEA(k, 1:2)); %Extract species properties
            isgas = str2double(CEA(k, 52));
            M = str2double(CEA(k, 53:65)); % (g/mol)
            h_form = str2double(CEA(k, 66:80)); % (J/mol)
            k=k+1;

            if numT == 0 % If data only at one T
                T = [];
                T(1) = str2double(CEA(k, 1:11));
                MAT.PRODUCTS.(name).T = T;
                numa = NaN;
                delh = NaN;
                k=k+1;
            elseif numT > 0 % If data extends over T range
                numa = str2double(CEA(k, 23));
                delh = str2double(CEA(k,66:80));
                n = [];
                a = [];
                T = [];
                for m = 1:numT % For each temperature range
                    T(m,1) = str2double(CEA(k, 1:11)); % (K)
                    T(m,2) = str2double(CEA(k, 12:22)); % (K)
                    for h = 1:numa % Store exponents
                        index = 19 + h*5;
                        n(m,h) = str2double(CEA(k, index:index+4));
                    end
                    k=k+1;
                    p = 0;
                    for h = 1:numa % Store coefficients
                        index = -15 + h*16;
                        if index > 80
                            index = index - 80;
                            if p == 0
                                k=k+1;
                                p = 1;
                            end
                        end
                        a_temp = CEA(k, index:index+15);
                        a_temp = strrep(a_temp, 'D', 'e');
                        a(m,h) = str2double(a_temp);
                    end
                    for h = 1:2 % Store integration constants
                        index = 33 + h*16;
                        b_temp = CEA(k, index:index+15);
                        b_temp = strrep(b_temp, 'D', 'e');
                        b(m,h) = str2double(b_temp);
                    end
                    k=k+1;
                end
                MAT.PRODUCTS.(name).T = T;
                MAT.PRODUCTS.(name).a = a;
                MAT.PRODUCTS.(name).n = n;
                MAT.PRODUCTS.(name).b = b;
            else
                error('thermo.inp is not fomatted properly. \n Problem at line %f', k);
            end
            MAT.PRODUCTS.(name).prop = [numT, isgas, M, h_form, numa, delh];
            k=k-1;
        end
    end  
    k=k+1;
end

% Loop through CEA looking for oxidizers
oxNames = fieldnames(MAT.OX);
k = 1;
while k < length(CEA)
    b1 = ~strcmp(CEA(k,1), {' ', '!', '#', '-'});
    b2 = ~strcmp(CEA(k,1:18), {'thermo            ', 'END PRODUCTS      ', 'END REACTANTS     '});
    
    if all(b1) && all(b2) % If the line isn't a comment, a space, or a keyword
        name = CEA(k,1:18); % Extract name of species
        name = strrep(name, '+', '_plus');
        name = strrep(name, '-', '_minus');
        name = matlab.lang.makeValidName(name);
%         if length(name) > 15
%             name = name(1:15); % Limit it to 15 characters because the CEA output is limited to 15 characters
%         end
        k=k+1;
        
        if any(strcmp(name, oxNames))
            numT = str2double(CEA(k, 1:2)); %Extract species properties
            M = str2double(CEA(k, 53:65)); % (g/mol)
            MAT.OX.(name).M = M;
            k=k+3*numT;
        end
    end  
    k=k+1;
end

% Loop through CEA looking for fuels
fuelNames = fieldnames(MAT.FUEL);
k = 1;
while k < length(CEA)
    b1 = ~strcmp(CEA(k,1), {' ', '!', '#', '-'});
    b2 = ~strcmp(CEA(k,1:18), {'thermo            ', 'END PRODUCTS      ', 'END REACTANTS     '});
    
    if all(b1) && all(b2) % If the line isn't a comment, a space, or a keyword
        name = CEA(k,1:18); % Extract name of species
        name = strrep(name, '+', '_plus');
        name = strrep(name, '-', '_minus');
        name = matlab.lang.makeValidName(name);
%         if length(name) > 15
%             name = name(1:15); % Limit it to 15 characters because the CEA output is limited to 15 characters
%         end
        k=k+1;
        
        if any(strcmp(name, fuelNames))
            numT = str2double(CEA(k, 1:2)); %Extract species properties
            M = str2double(CEA(k, 53:65)); % (g/mol)
            MAT.FUEL.(name).M = M;
            k=k+3*numT;
        end
    end  
    k=k+1;
end

% 
% % Generate thermodynamic property lookup table for products
% for k = 1:length(prodNames) % Iterate through each product species
%     tp = thermo.(prodNames{k}); % Set up temporary variable to improve legibility
%     M = tp.prop(3)/1000; % (kg/mol)
%     Rk = 8.3144598/M; % (J/kg*K)
%     
%     T = tp.T(1):1:tp.T(end); % Change the middle number to define the fineness of the temperature step
%     
%     ii = 1;
%     mm = zeros(size(T));
%     for kk = 1:length(T) % Sets mm (index for temperature range) to correct range for range T
%         if T(kk) <= tp.T(ii,2)
%             mm(kk) = ii;
%         else
%             ii = ii + 1;
%             mm(kk) = ii;
%         end
%     end
%     
%     cp = zeros(size(T));
%     h = zeros(size(T));
%     for kk = 1:length(T)
%         cp(kk) = Rk*sum( tp.a(mm(kk),:) .* (T(kk) .^ (tp.n(mm(kk),:))) ); % (J/kg*K)
%         h_temp = (tp.a(mm(kk),:)./(tp.n(mm(kk),:)+1)) .* (T(kk) .^ (tp.n(mm(kk),:)+1)); % Uses the formula for integration of a polynomial, since delh = integral(cp)dT
%         h_temp(~isfinite(h_temp)) = tp.a(mm(kk), ~isfinite(h_temp)).* log(T(kk)); % Replaces the infinite term with ln(T) since polynomial integration formula doesn't work for 1/T
%         h(kk) = Rk*(sum(h_temp) + tp.b(mm(kk),1)); % (J/kg)
%     end
%     
%     MAT.PRODUCTS.(prodNames{k}).T = T;
%     MAT.PRODUCTS.(prodNames{k}).cp = cp;
%     MAT.PRODUCTS.(prodNames{k}).h = h;
%     MAT.PRODUCTS.(prodNames{k}).M = tp.prop(3);
% end
% 
% % Add molecular weight to oxidizer
% oxNames = fieldnames(MAT.OX);
% for k = 1:length(oxNames)
%     MAT.OX.(oxNames{k}).M = thermo.(oxNames{k}).prop(3);
% end
% 
% % Add molecular weight to fuel
% 
% for k = 1:length(oxNames)
%     MAT.FUEL.(fuelNames{k}).M = thermo.(fuelNames{k}).prop(3);
% end

% Save as '.mat' file
save([folder file '.mat'],'-struct','MAT')
fclose(fileID);
end

