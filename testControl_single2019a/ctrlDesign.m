function Controller = ctrlDesign(patient, time_vec, Food)
% Use this template to design an close-loop system controller to stablize the 
% patient sugar response

%% system identification (open loop sim)
% The input response is loaded here and used to simulate the patient to produce 
% the step response. Feel free to alter this section as needed to try different
% types of inputs that may help with the identification process
[TF,IC] = sysID(patient); % REPLACE THIS FUNCTION (AT THE BOTTOM) WITH YOURS

%% system identification (close loop sim)
% Simulate the open loop response of the generated patient
Controller = tf(0); % setting a null controller to get the closed loop response without a controller
Sugar_closeloop = closedLoopSim(patient,Food,Controller);

% Get Sugar values at time_vec time. This is basic linear interpolation and
% is nessesary because Simulink does not guarantee Sugar.Time will equal time_vec
sugar_vec_closeloop = interp1(Sugar_closeloop.Time,Sugar_closeloop.Data,time_vec,'linear');

%% controller design
s = tf('s');
Controller = (-0/(s+2));
end

function [TF, IC] = sysID(patient)
% Use this template to design an open-loop system identification routine given
% the step time response of the patient. 

%% input response
% The input response is loaded here and used to simulate the patient to produce 
% the step response. Feel free to alter this section as needed to try different
% types of inputs that may help with the identification process
[time_vec, Food, InsulinRate] = inputVector();

% Simulate the open loop response of the generated patient
Sugar = openLoopSim(patient,Food,InsulinRate);

% Get Sugar values at time_vec time. This is basic linear interpolation and
% is nessesary because Simulink does not guarantee Sugar.Time will equal time_vec
sugar_vec = interp1(Sugar.Time,Sugar.Data,time_vec,'linear');

%% system identification
%Find the peaks in the data, this is when the slope changes sign
[maxPKS,maxLOCS] = findpeaks(sugar_vec,time_vec);
[minPKS,minLOCS] = findpeaks(-sugar_vec,time_vec);
minPKS = -minPKS; 

%Find the max 
max_val = max(sugar_vec);
%Find steady state value
steady = sugar_vec(end);
s = tf('s');

if length(minPKS)<1
    damp = 7.5;
    wn = 2;
    a = 4/570;
else
    
    % OS
    OS = (steady-minPKS(1))/steady;
    %time of first peak (in our case a local minima)
    tp = minLOCS(1);
        damp = abs(log(OS)/sqrt( pi^2 + log(OS)^2 )) * 0.8;
        wn = (pi/tp) / (sqrt(1-damp^2))*0.85;    
        a = wn*damp*steady*0.33;

    if OS < 0 || minLOCS(1)>600
        damp = 7.5;
        wn = 2;
        a = 4/570;
    elseif length(maxPKS)>1
        if maxPKS(2)>(steady+1)
        damp = 0.85;
        wn = 0.75*(pi/tp) / (sqrt(1-damp^2));
        a = wn*damp*steady;
        end
    end
end
    %Transfer function that is used for all 'types' of systems
    TF = -(max_val - steady)*a*(wn^2 + 10^(-10)) / ((s+a)*( (s^2) + (2*wn*damp*s) + (wn^2))); 
    %Produce initial condition (offset from zero)
    IC = max_val; %changed the start value to the max value

end
%Useless Code
% 
% val90 = (steady) +(max_val - steady)*0.1;
% val10 = (steady) +(max_val - steady)*0.9;
% array90 = logical(sugar_vec<val90);
% array10 = logical(sugar_vec<val10);
% t10 = find(array10, 1, 'first');
% t90 = find(array90, 1, 'first');
% tr = (t90-t10)*1.6;
% 
%ts = +-2%
% ts = 0;
% tsarray = logical((steady-0.02*steady)<sugar_vec & sugar_vec<(steady+0.02*steady));
% for i = 0:length(tsarray)-1
%     i_rev = length(tsarray)-i;
%     if(tsarray(i_rev) == 0)
%         if (i_rev>ts)
%             ts = i_rev;
%         end
%     end
% end