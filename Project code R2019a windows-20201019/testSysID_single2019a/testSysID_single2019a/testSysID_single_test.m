clearvars;
close all;
clc;
marks = 0;
total_patients = 10; %Run 10 patients
Total_mark = 0;

for i = 1:total_patients
% prepare input signal
[time_vec, Food, InsulinRate] = inputVector();
insulin_step = ones(size(time_vec))';

% Generate a new random patient and simulate the open loop response of the
% generated patient
patient = genPatient();

% Simulate the actual patient, then interpolate since Simulink does not
% guarantee Sugar.Time will equal time_vec
Sugar = openLoopSim(patient,Food,InsulinRate);
sugar_vec = interp1(Sugar.Time,Sugar.Data,time_vec,'linear');

% Simulate the open loop response of the system id process
[TF,IC] = sysID(patient);
Y = step(TF,time_vec);
id_resp = Y+IC;

% Simulate the reference system
[TF_ref, IC_ref] = referenceID(patient);
Y_ref = step(TF_ref,time_vec);
ref_resp = Y_ref+IC_ref;

% Preparing the data for plotting
time = time_vec/60; % convert to hours
patient_sugar_resp = sugar_vec(:);
id_sugar_resp = id_resp(:);
ref_sugar_resp = ref_resp(:);

% Set up anon function for RMSE calculation
rmseFct = @(x, y) sqrt(sum((normVector(x - y)).^2)/(size(x, 1)));
rmse_id = rmseFct(patient_sugar_resp, id_sugar_resp);
rmse_ref = rmseFct(patient_sugar_resp, ref_sugar_resp);

%Mark Scheme 
mark = 10 - 10*(rmse_id /rmse_ref);
marks = marks + mark;
Order_TF = max(length(TF.Numerator{1,1})-1,length(TF.Denominator{1,1})-1);
Scaling = 5/max(5,Order_TF);
Patient_mark = min(10,Scaling*mark);
Total_mark = Total_mark + Patient_mark;
fprintf("Patient No. = %d\n Mark (Before scaling) = %f\n Mark (After scaling) = %f\n\n", i,mark,Patient_mark);

%plotSysId(time, patient_sugar_resp, ref_sugar_resp, id_sugar_resp, rmse_id, rmse_ref);
end
fprintf("Total mark (out of 100) = %d\n\n",Total_mark);
