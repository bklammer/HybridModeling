function [MAT] = CEAtoMATLAB3(file, folder)
%CEAtoMATLAB takes a CEA output file and converts it to a .mat file
% 'folder' is a string specifying a folder name to look for file in
% 'file' is a desired file name for the saved .mat file

% Handle different inputs
if nargin == 0
    folder = '';
    file = 'MuleSim3CEA';
elseif nargin == 1
    folder = '';
end

%% READING COMBUSTION PROPERTIES FILE %%
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
k = 1; % Counter for length of array CEA
m = 1; % Index for different OF ratios
n = 1; % Index for different pressures
O = 1; % Check to see if oxidizer data has already been recorded
F = 1; % Check to see if fuel data has already been recorded
while k < length(CEA)

    if O == 1 && strcmp(CEA{k}, 'OXIDANT') % Extract oxidizer data, only once
        k=k+1;
        while ~strcmp(CEA{k}, 'FUEL') && k < length(CEA)
            if strcmp(CEA{k}, 'OXIDANT')
                k=k+1; % If another oxidizer detected, increment
            end
            oxName = CEA{k}; % Record name of oxidizer and convert it to a valid MATLAB variable name
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
        O = 0; % Once oxidizer data extracted, never use this loop again
    end
    
    if F == 1 && strcmp(CEA{k}, 'FUEL') % Extract fuel data, only once
        k=k+1;
        while ~strcmp(CEA{k}, 'O/F=') && k < length(CEA)
            if strcmp(CEA{k}, 'FUEL') 
                k=k+1; % If another fuel is detected, increment
            end
            fuelName = CEA{k}; % Record name of fuel and convert it to a valid MATLAB variable name
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
        F = 0; % Once fuels extracted, don't look for more fuels
    end
    
    if strcmp(CEA{k}, 'O/F=') % Store O/F data
        k=k+1;
        OF = str2double(CEA{k});
        if MAT.OF(m) ~= OF && MAT.OF(m) ~= 0  
            m=m+1; % Increment OF index if new OF ratio detected
            MAT.OF(m) = OF;
            n = 1; % Reset pressure index n if new OF ratio detected
        else
            MAT.OF(m) = OF; % If first OF ratio is detected, do not increment m, do not reset n
        end
    end
    
    if m == 1 && strcmp(CEA{k}, 'BAR') % Store pressure data only for the first OF ratio
        k=k+1;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.p(n_temp) = str2double(CEA{k})*100000; % Convert pressure to Pa and record
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'T,') % Store temperature data
        k=k+2;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.T(m, n_temp) = str2double(CEA{k}); % Record temperature
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'RHO,') % Store density data
        k=k+3;
        n_temp = n;
        while k < length(CEA) && ~strcmp(CEA{k}, 'H,')
            if isnan(str2double(CEA{k})) % if the exponent is negative (i.e. '2.7181-1')
                rho_temp = split(CEA{k},'-');
                rho = str2double(rho_temp{1});
                expo = -str2double(rho_temp{2});
            else
                rho = str2double(CEA{k});
                k=k+1;
                expo = str2double(CEA{k});
            end
            MAT.rho(m, n_temp) = rho*10^expo; % Record density
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'H,') % Store enthalpy data
        k=k+2;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.h(m, n_temp) = str2double(CEA{k})*1000; % Convert to J/kg and record enthalpy
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'U,') % Store internal energy data
        k=k+2;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.u(m, n_temp) = str2double(CEA{k})*1000; % Convert to J/kg and record internal energy
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'M,') % Store molar mass data
        k=k+2;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.M(m, n_temp) = str2double(CEA{k})/1000; % Convert to kg/mol and record molar mass
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'Cp,') % Extract specific heat at constant pressure data
        k=k+2;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.cp(m, n_temp) = str2double(CEA{k})*1000; % Convert to J/kg*K and record specific heat
            n_temp=n_temp+1;
            k=k+1;
        end
    end
    
    if strcmp(CEA{k}, 'GAMMAs') % Extract specific heat ratio
        k=k+1;
        n_temp = n;
        while k < length(CEA) && ~isnan(str2double(CEA{k}))
            MAT.k(m, n_temp) = str2double(CEA{k}); % Record specific heat ratio
            n_temp=n_temp+1;
            k=k+1;
        end
        n = n_temp; % Increment n since there might be more values at the same OF ratio lower in the file
    end
    
    k=k+1;
end

%% READING THERMODYNAMIC DATABASE FILE %%
% Import 'thermo.inp' file to char array
[fileName, filePath] = uigetfile([folder, '*.inp'], 'Please select a CEA thermodynamic database');
if isequal(fileName, 0)
    error('User selected cancel')
end
fullFileName = [filePath, fileName];
C = regexp(fileread(fullFileName), '\r?\n', 'split');
CEA = char(C{:});
CEA = CEA(:,1:80);

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

% Save as '.mat' file
[file,path] = uiputfile([file, '.mat']); % User input file name
if file == 0
    error('User selected cancel')
end
filename = fullfile(path, file);
save(filename,'-struct','MAT')
fclose(fileID);
end

