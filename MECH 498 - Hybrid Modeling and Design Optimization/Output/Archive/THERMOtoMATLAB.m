function [MAT] = THERMOtoMATLAB(folder, file)
%THERMOtoMATLAB takes a CEA thermo input file and converts it to a .mat file
% 'folder' is a string specifying a folder name to look for file in
% 'file' is a desired file name for the saved .mat file

% Handle different inputs
if nargin == 0
    folder = '';
    file = 'MuleSimTHERMO';
elseif nargin == 1
    file = 'MuleSimTHERMO';
end

% Import 'thermo.inp' file to char array
[fileName, filePath] = uigetfile([folder, '*.inp']);
if isequal(fileName, 0)
    error('User selected cancel')
end
fullFileName = [filePath, fileName];
C = regexp(fileread(fullFileName), '\r?\n', 'split');
CEA = char(C{:});
CEA = CEA(:,1:80);

MAT = struct();
k = 1;
while k < length(CEA)
    b1 = ~strcmp(CEA(k,1), {' ', '!', '#', '-'});
    b2 = ~strcmp(CEA(k,1:18), {'thermo            ', 'END PRODUCTS      ', 'END REACTANTS     '});
    
    if all(b1) && all(b2) % If the line isn't a comment, a space, or a keyword
        name = CEA(k,1:18); % Extract name of species
        name = strrep(name, '+', '_plus');
        name = strrep(name, '-', '_minus');
        name = matlab.lang.makeValidName(name);
        if length(name) > 15
            name = name(1:15); % Limit it to 15 characters because the CEA output is limited to 15 characters
        end
        k=k+1;
        
        numT = str2double(CEA(k, 1:2)); %Extract species properties
        isgas = str2double(CEA(k, 52));
        M = str2double(CEA(k, 53:65)); % (g/mol)
        h_form = str2double(CEA(k, 66:80)); % (J/mol)
        k=k+1;
            
        if numT == 0 % If data only at one T
            T = [];
            T(1) = str2double(CEA(k, 1:11));
            MAT.(name).T = T;
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
            MAT.(name).T = T;
            MAT.(name).a = a;
            MAT.(name).n = n;
            MAT.(name).b = b;
            
        else
            error('thermo.inp is not fomatted properly. \n Problem at line %f', k);
        end
        
        MAT.(name).prop = [numT, isgas, M, h_form, numa, delh];
        k=k-1;
    end  
    k=k+1;
end

% Save as '.mat' file
save([folder file '.mat'],'MAT')

end

