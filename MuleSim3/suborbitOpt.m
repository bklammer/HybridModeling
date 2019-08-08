function [p] = suborbitOpt(p)
%-------------------------------------------------------------------------%
% suborbit function
% Michael Pearson - 2016
%
% Modified by Benjamin Klammer for design optimization of Ramses-1 - 2019
%-------------------------------------------------------------------------%
% Function Description
%-------------------------------------------------------------------------%
% Calculates the altitude, velocity and acceleration of a rocket
%
% Input is struct 'p', which must contain fields:
%   'dia' containing the outer diameter of the rocket in m
%   'mass' containing the total mass of the rocket in kg
%   'launchAlt' containing the launch elevation in m
%   'p_amb' containing the ambient pressure used in MuleSim3 in Pa
%   'A_th' containing the nozzle throat area in m^2
%   'AR' containing the nozzle area ratio
%   'TC' containing an array with burn time in column 1, thrust data in
%    column 2, mass in column 3, and with the first row containing total
%    impulse in column 1, average thrust in column 2, and total mass
%    including propellant and case in column 3
%-------------------------------------------------------------------------%
%                                  _
%                                 / \
%                                 | |
%                                 | |
%                                 | |
%                                /| |\
%                               / | | \
%                               |/| |\|
%                               ' *** '
%                                  *
%-------------------------------------------------------------------------%


CdA = p.CdA;
dia = p.dia;

area = pi*(1.05*dia/2)^2;       % factor of 1.05 accounts for fins
time = 0;           % time accounting
dt = 0;             % time delta
%gravloss = 0;       % determines total drag loss to delta-v

onpower = 0;        % for plot to show when engine is firing
stagesep = 0;       % for plot to show when rocket is staging

launchAlt = p.launchAlt; % Launch altitude in meters - must be greater than 0
p_amb = p.p_amb; % Ambient pressure used in MuleSim3 for thrust correction
A_exit = p.AR*p.A_th; % Nozzle exit area for thrust correction

TC = p.TC; % Thrust/mass curve
mass = p.mass; % Total mass of rocket
[alt1,vel1,accel1,M1] = power(launchAlt,0,0,TC,CdA);
[alt2,vel2,accel2,M2] = coast(alt1(end),vel1(end),accel1(end),CdA);

alt = [alt1;alt2]; vel = [vel1;vel2]; accel = [accel1;accel2]; M = [M1;M2];% Concatenate power and coasting arrays
alt = alt-launchAlt; % Subtract launch altitude to get AGL altitude
p.alt = alt;
p.vel = vel;
p.accel = accel;
p.M = M;
p.apogee = alt(end);

%-------------------------------------------------------------------------%
% Post processing for data output
%-------------------------------------------------------------------------%
% onpower(onpower == 0) = NaN; onpower(length(alt),1) = NaN;
% stagesep(stagesep == 0) = NaN; stagesep(length(alt),1) = NaN;
% SOoutput(alt,vel,accel,onpower,stagesep,time);


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                               /  |
% Sub Functions         ,------'   '-------..
%                      (____________________-+
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Rocket motor producing thrust
%-------------------------------------------------------------------------%
    function [alt,vel,accel,M] = power(alt,vel,accel,TCurve,CdA)
        dt = TCurve(end,1)/length(TCurve); % Assumes evenly spaced thrust curve
        plen = length(onpower);
        
        for i=1:(length(TCurve)-2)
            %%% calculate atmospheric conditions %%%
            [~,a,p_act,rho] = atmosisa(alt(i));
            
            %%% calculate drag force %%%
            M(i,1) = vel(i)/a; % Mach number
            [~,index] = min(abs(CdA(:,1)-M(i,1))); % Find the index of the closest mach number in Cd table to actual mach number
            drag = CdA(index,2)*area*0.5*rho*vel(i)^2;
            
            %%% calculate gravitational acceleration %%%
            g = grav(alt(i)); % gravity decreases as you get farther from earth
            
            %%% correction for sea-level thrust curves %%%
            Thrust(i) = TCurve(i+1,2) + (p_amb - p_act)*A_exit;
            
            %%% calculate acceleration, velocity and altitude %%%
            accel(i+1,1) = (Thrust(i)-drag)/mass-g;
            vel(i+1,1) = accel(i)*dt+vel(i);
            alt(i+1,1) = vel(i)*dt+((accel(i)*dt^2)/2) +alt(i);
            
            if alt(i+1,1) < launchAlt   % check if hit ground
                alt(i+1,1) = launchAlt;
                vel(i+1,1) = 0;
                accel(i,1) = 0;
            end
            
            onpower(i+plen+1,1) = alt(i+1,1);       % for plotting
            stagesep(i+plen+1,1) = NaN;             % for plotting
            
            mass = mass - (TCurve(i+1,3)-TCurve(i+2,3)); % Michael's motor data mass data is in grams
            time = time + dt;
 %           gravloss = gravloss + (drag*dt/mass);
        end
        
        index = cast(1+p.burn, 'int16'); % Turn burn time float into integer
        m_out = TCurve(2,3) - TCurve(index,3);
        p.I_sp = trapz(Thrust(1:(index-1)))*dt/(9.81*m_out); % (s)
    end % end power

%-------------------------------------------------------------------------%
% Rocket Coasting to apogee
%-------------------------------------------------------------------------%
    function [alt,vel,accel,M] = coast(alt,vel,accel,CdA)
        i = 1;
        plen = length(onpower);
        
        while vel(i) > 0
            %%% calculate atmospheric conditions %%%
            [~,a,~,rho] = atmosisa(alt(i));
            
            %%% calculate drag force %%%
            M(i,1) = vel(i)/a; % Mach number
            [~,index] = min(abs(CdA(:,1)-M(i,1)));
            drag = CdA(index,3)*area*0.5*rho*vel(i)^2;
            
            %%% calculate gravitational acceleration %%%
            g = grav(alt(i)); % must take into account
            
            %%% calculate acceleration, velocity and altitude %%%
            accel(i+1,1) = -drag/mass-g;
            vel(i+1,1) = accel(i)*dt+vel(i);
            alt(i+1,1) = vel(i)*dt+((accel(i)*dt^2)/2) +alt(i);
            
            if alt(i+1,1) < launchAlt   % check if hit ground
                alt(i+1,1) = launchAlt;
                vel(i+1,1) = 0;
                accel(i,1) = 0;
            end
            
            onpower(i+plen,1) = NaN;          % for plotting
            stagesep(i+plen,1) = NaN;         % for plotting
            
            i = i + 1;
            time = time + dt;
 %           gravloss = gravloss + (drag*dt/mass);
        end
    end % end coast

end % end main program
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%% Calculates gravitational acceleration at altitude
function g = grav(alt)
g = (6.67408e-11*5.972e24)/(6.371e6+alt)^2;
end

%% Outputs data and creates plots
function SOoutput(alt,vel,accel,onpower,stagesep,time)
data(1,:) = cellstr(['m    ';'m/s  ';'m/s^2';'s    ']);
data(2,:) = num2cell([max(alt),max(vel),max(accel),time]);

% *** To output performance table, uncomment this section ***
array2table(data,'VariableNames'...
        ,{'Apogee' 'Vmax' 'Accelmax' 'TimeToApogee'},'RowNames'...
        ,{'Units','Value'})

t = linspace(0,time,length(alt));
plot(t,vel,t,accel);
xlim([0,t(end)]);
ylim([-500,max(alt)*1.1]);
xlabel('Time (s)');
ylabel('Altitude (m)');
hold on
% comet(t,alt)
plot(t,alt)
plot([0,t(end)],[100000,100000],':')
plot(t,onpower,'Linewidth',2)
plot(t,stagesep,'c','Linewidth',2)
hold off

end