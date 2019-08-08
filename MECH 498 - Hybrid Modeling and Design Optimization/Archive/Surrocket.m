function [out, p] = Surrocket(p)
%SURROCKET Models a hybrid rocket from design to apogee
%   Detailed explanation goes here


%% INPUTS %%
% Load input files
input = readtable('MuleSim3INPUT.xlsx');
CEA = load('MuleSim3CEA.mat');
N2Osat = load('N2Osat.mat');

% Pre-declare variables before they're poofed into existence by eval
[del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, d_f, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    zeta_d, zeta_cstar, zeta_CF, d_th, p_atm, A_ratio_nozzle, graph, save] = deal(NaN);

for k = 1:length(input.Symbol)
    if ~exist(input.Symbol{k}, 'var')
        error('Variable ''%s'' in the input spreadsheet is not declared in the MATLAB code', input.Symbol{k})
    end
    eval([input.Symbol{k} '= input.Value(' num2str(k) ');']); % Read input data from spreadsheet
end

tic; % Start recording the simulation runtime

design = [del_time, time_max, V_tank, p_tank_init, m_ox_tank_init, n_inj, d_inj, ...
    Cd_inj, d_f, rho_f, a, n, L, d_port_init, T_cc_init, p_cc_init, ...
    zeta_d, zeta_cstar, zeta_CF, d_th, p_atm, A_ratio_nozzle];

p.design = design;
p.CEA = CEA;
p.N2Osat = N2Osat;

%% MODEL CALCS

% Hybrid Motor Model
p = MuleSim3Opt(p);

% Structural Mass and Diameter Estimation
p.dia = p.d_port_f + 0.1;
p.mass_structural = 0.434866*p.mass_prop + 259.0147*p.dia - 13.2621;
p.mass = p.mass_structural + p.mass_prop;

% Trajectory Model
p = suborbitOpt(p);

% Output
out = p.apogee;

toc;

end
                        
% 
%                                         @@&              
%                                  %%   @  @@  @           
%                              @      *@     @             
%                    @@@@%   @      @@     @/#@@           
%                  .@@@     @ /@  @,            @          
%                  .@@@     (   @              %@          
%                   @@@@    @       @@   @    %@           
%                     @@@@ .%,@.  @ *   @    @             
%                       @@@@( @ @@ @   /  @/   @@          
%                      @  *@@@@@     #         @@/         
%                     @       @@@@@@@@.     ,@@@@          
%                     @     @      %@@@@@@@@@@@.           
%                    & @(                                  
% 



%                                             &@&.                       
%                                          @@      /@                    
%                               %@@@@@@,  @&    @%   %(                  
%                           (@%         @@@        @                     
%              ,&@@@@@@@@@@@.         @@&         @#                     
%          *@@@@@@&      @/         @@,       ,&,  /@@@.                 
%         @@@@@%        @    &@@@@@@.                 @@%                
%        #@@@@@        @..@*    @@                     @@                
%        *@@@@@        @,    (@/                      &@,                
%         @@@@@@          @@.         *@@@@@,        #@#                 
%          @@@@@@    (@#           #@@      @       @@.                  
%            @@@@@@  .&@@@@@@    @@ @      @/     /@&                    
%             #@@@@@@.    #@   &@  @      @     @@/  #@,                 
%               .@@@@@@@. @@  @@@  @    @.   @@%     @@@%                
%               @  @@@@@@@@@ % @  ,   @%@@@*         #@@@                
%             /#      %@@@@@@@@@.                    @@@@/                       
%            /%           @@@@@@@@@@@@,           (@@@@@@                
%             @          *@.  *@@@@@@@@@@@@@@@@@@@@@@@@@                 
%            @/      .@@            ,&@@@@@@@@@@@@@@@                    
%           @    @@,                                                     
%          @@%                                                           
