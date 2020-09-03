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

% Here are some potentially useful functions:
% - findpeak
% - min/max

%Find the peaks in the data, this is when the slope changes sign
[maxPKS,maxLOCS] = findpeaks(sugar_vec,time_vec);
[minPKS,minLOCS] = findpeaks(-sugar_vec,time_vec);
minPKS = -minPKS; 
Ts = max([max(maxLOCS), max(minLOCS)]); %found the last local minima/maxima

%Find the min and max 
min_val = min(sugar_vec);
max_val = max(sugar_vec);

%In the example we ignore the peaks, min and max and just output a first order 
%response shifted by 160 
%Produce first order transfer function as output
s = tf('s');
% TF = -(max_val-min_val)*(4/(600*s+4)); %this one is like the reference
% ^^ just changed the value for steady state and simplified the transfer
% func

%need to adjust poles...
% TF = -(max_val-min_val)*(4/(minLOCS(1)*s+4))%this one drops too early
% TF = -(max_val-min_val)*(4/((600*s+4)*(600*s-4))); %not good haha
wn = 1.8/(minLOCS(1)/2.5); %~0.0054
damp = 0.9;
TF = -(max_val-min_val)*((wn^2)/( (s^2) + (2*wn*damp*s) + (wn^2))); 
%^^ works okay for some things (better than ref) 
% except when second oscillation is quite large

%Produce initial condition (offset from zero)
IC = max_val; %changed the start value to the max value

end