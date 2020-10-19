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
Controller = tf(-0/(s+1));
end

function [TF, IC] = sysID(patient) % update this function as appropriate
%% input response
[time_vec, Food, InsulinRate] = inputVector();
Sugar = openLoopSim(patient,Food,InsulinRate);
sugar_vec = interp1(Sugar.Time,Sugar.Data,time_vec,'linear');

%% system identification
[PKS,LOCS] = findpeaks(sugar_vec,time_vec);
min_val = min(sugar_vec);
max_val = max(sugar_vec);

s = tf('s');
TF = -80*(4/(10*60))/(s+4/(10*60));
IC = 160;
end