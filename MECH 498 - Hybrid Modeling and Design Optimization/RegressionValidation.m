function [d_port_final, OF_average] = RegressionValidation()
%REGRESSIONVALIDATION outputs error for given motors

% Benjamin Klammer - 2019
% All calculations performed in metric units

%% INPUTS %%
del_time = 0.01; % (s) time step
time_max = [8.3, 8.25, 8.15, 8.4]; % (s) burn time
m_dot_ox_in = [4.44, 4.43, 4.42, 4.43]; % (kg/s)
a = [9.36e-5, 9.24e-5, 9.10e-5, 9.36e-5]; % Average of four Stanford/Nasa tests, gives r_dot in m/s
n = 0.62;
L = [1.149, 1.149, 1.148, 1.148]; % (m)
d_port_init = [0.0893, 0.1001, 0.1030, 0.1138]; % (m)
rho_f = 900; % (kg/m^3) Density of fuel

graph = 0;

tic; % Start recording the simulation runtime


%% MODEL CALCULATIONS %%

for p = 1:length(time_max) % Run all four motors back to back
    t = 1;
    r_cc = d_port_init(p)/2;
    while 1 % Run until break
        A_cc(t) = pi()*r_cc(t)^2; % Port geometry
        G(t) = m_dot_ox_in(p)/A_cc(t);
        del_r_cc(t) = a(p)*G(t)^n * del_time; % Regression rate law
        m_dot_f(t) = 2*pi()*r_cc(t)*L(p) * rho_f * (del_r_cc(t)/del_time); % Fuel mass flow rate from geometry
        
        m_dot_f_temp = 0;
        k = 0;
        while abs(m_dot_f_temp - m_dot_f(t)) > 0.01*m_dot_f(t) && k < 100 % until m_dot_f converges
            m_dot_f_temp = m_dot_f(t);
            m_dot_cc(t) = m_dot_f(t) + m_dot_ox_in(p); % Total mass flow in
            G(t) = (m_dot_ox_in(p) + m_dot_cc(t))/(2*A_cc(t)); % Average mass flux accross entire port
            del_r_cc(t) = a(p)*G(t)^n * del_time; % Regression rate law using updated G
            m_dot_f(t) = 2*pi()*r_cc(t)*L(p) * rho_f * (del_r_cc(t)/del_time); % Fuel mass flow rate from geometry
            k = k+1;
        end

        OF(t) = m_dot_ox_in(p)./m_dot_f(t); % Calculate Oxidizer-Fuel Ratio

        if t < time_max(p)/del_time % If less than burn time
            r_cc(t+1) = r_cc(t) + del_r_cc(t);
            t = t+1;
        else
            break % Exit loop if no more time or oxidizer
        end
    end
    d_port_final(p) = 2*r_cc(end); % Save final port diameter
    OF_average(p) = mean(OF);
end

time = 0:del_time:(del_time*(t-1)); % Create a time vector
time_sim = toc; % Save how long the simulation took (minus loading data and plotting figures)


%% PLOTS %%
if graph == 1 % Graph if indicated
    figure
    hold on
    plot(time, r_cc.*1000) % (mm)
    title('Stanford Fuel Grain Radius')
    xlabel('Time (s)')
    ylabel('Fuel Grain Radius (mm)')
    axis([0 9 0 90])
    hold off
end


end

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

