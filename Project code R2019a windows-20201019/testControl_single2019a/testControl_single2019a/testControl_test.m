clear all;
close all;
clc;
steadystate_desired(1) = 3.5*18; % convert mmol/l to mg/dl 
steadystate_desired(2) = 7*18; % convert mmol/l to mg/dl 
peak_dangerous(1) = 2.2*18; % convert mmol/l to mg/dl 
peak_dangerous(2) = 16.6*18; % convert mmol/l to mg/dl 
marks = 0;
Total_mark = 0;
total_patients = 10;


for i = 1:total_patients %Run multiple patients at once
tic;
%Generate a patient
patient = genPatient();

% comment whichever appropriate one as needed 
%[time_vec, food] = foodVector_fasting(); % simulate fasting response
[time_vec, food] = foodVector_3meals(); % simulate 3 meals

%Create a Controller
Controller = ctrlDesign(patient, time_vec, food);

%Simulate closed loop system
Sugar = closedLoopSim(patient,food,Controller);

time = Sugar.Time/60;
patient_sugar_resp = Sugar.Data(:);

%Plot results
fig = plotCtrlDesign(time, patient_sugar_resp, steadystate_desired, peak_dangerous);
Run_time = toc;

% Mark Scheme
    mark = 10*(sum(logical(patient_sugar_resp > steadystate_desired(1)) & logical(patient_sugar_resp < steadystate_desired(2)))/length(patient_sugar_resp));
     
    Lower_threshold = sum(logical(patient_sugar_resp < 40)); %Checking if any value is less than 40
    Upper_threshold = sum(logical(patient_sugar_resp > 290)); %Checking if any value is greater than 290
    
    if sum(Lower_threshold) > 0 % if any sample is less than 40mg/dl, mark is 0 for that patient
       mark = 0; 
    elseif sum(Upper_threshold) > 0 % if any sample is greater than 290mg/dl, mark is 0 for that patient
       mark = 0;
    elseif Run_time > 60  %if run time is greater than a minute, mark is 0 for that patient
       mark = 0;
    end
    
    marks = marks + mark;
    Order_TF = max(length(Controller.Numerator{1,1})-1,length(Controller.Denominator{1,1})-1);
    Scaling = 5/max(5,Order_TF);
    Patient_mark = min(10,Scaling*mark);
    Total_mark = Total_mark + Patient_mark;
    fprintf("Patient No. = %d\n Mark (Before scaling) = %f\n Mark (After scaling) = %f\n\n", i,mark,Patient_mark);
end
fprintf("Total marks (out of 100) = %d\n\n",Total_mark);