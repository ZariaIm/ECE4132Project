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
[PKS,LOCS] = findpeaks(sugar_vec,time_vec);

%Find the min and max 
min_val = min(sugar_vec);
max_val = max(sugar_vec);

%In the example we ignore the peaks, min and max and just output a first order 
%response shifted by 160 
%Produce first order transfer function as output
s = tf('s');
TF = -(max_val-min_val)*(4/(10*60))/(s+4/(10*60));

%Produce initial condition (offset from zero)
IC = max_val;

%notes: steady state value = min, initial value/offset = max, 
%atm my first order system is as accurate as the reference but basically
%the same as it.
%A second order system should make it more accurate.

end