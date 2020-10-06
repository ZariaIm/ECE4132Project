close all; clc;
steadystate_desired(1) = 3.5*18; % convert mmol/l to mg/dl 
steadystate_desired(2) = 7*18; % convert mmol/l to mg/dl 
peak_dangerous(1) = 2.2*18; % convert mmol/l to mg/dl 
peak_dangerous(2) = 16.6*18; % convert mmol/l to mg/dl 

%Generate a patient
patient = genPatient();

% comment whichever appropriate one as needed 
[time_vec1, food1] = foodVector_fasting(); % simulate fasting response
[time_vec2, food2] = foodVector_3meals(); % simulate 3 meals

%Create a Controller
Controller1 = ctrlDesign(patient, time_vec1, food1);
Controller2 = ctrlDesign(patient, time_vec2, food2);

%Simulate closed loop system
Sugar1 = closedLoopSim(patient,food1,Controller1);
Sugar2 = closedLoopSim(patient,food2,Controller2);

time1 = Sugar1.Time/60;
time2 = Sugar2.Time/60;
patient_sugar_resp1 = Sugar1.Data(:);
patient_sugar_resp2 = Sugar2.Data(:);

%Plot results
fig1 = plotCtrlDesign(time1, patient_sugar_resp1, steadystate_desired, peak_dangerous);
xlim([0 24])
title('Fasting')

fig2 = plotCtrlDesign(time2, patient_sugar_resp2, steadystate_desired, peak_dangerous);
xlim([0 24])
title('3meals')