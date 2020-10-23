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
p1 = -5;
p2 = -0.5;
z1 = -1.2;
Controller = tf( (-(s-z1)) / ((s-p1)*(s-p2)) );
end

function [TF, IC] = sysID(patient) % update this function as appropriate
%% input response
[time_vec, Food, InsulinRate] = inputVector();
Sugar = openLoopSim(patient,Food,InsulinRate);
sugar_vec = interp1(Sugar.Time,Sugar.Data,time_vec,'linear');

%% system identification
% Produce initial condition (offset from zero)
IC = sugar_vec(1); %changed the start value to the max value
% simulate the reference as comparison
s = tf('s');
min_val = min(sugar_vec);
max_val = max(sugar_vec);
TF_ref = -(max_val-min_val)*(4/(600*s+4));
Y_ref = step(TF_ref,time_vec);
ref_resp = Y_ref+IC;
ref = ref_resp(:);
sugar = sugar_vec(:);

% calculate the rmse of the reference
rmseFct = @(x, y) sqrt(sum((normVector(x - y)).^2)/(size(x, 1)));
rmse_ref = rmseFct(sugar, ref);
% make function handles for the id
num = @(a) [10*a(1) -0.5*a(2) + 0.0001*a(3) (6e-6)*a(4) (1e-5)*a(5) (5e-9)*a(6) 0];
den =@(a) [1, (0.0001)*a(7), (5e-7)*a(9), (-3e-7)*a(8) 0 1.782e-15];
TF_temp =@(a) tf(num(a), den(a));
id =@(a) (step(TF_temp(a) ,time_vec) + IC);
rmse_id = @(a) rmseFct(sugar, id(a));
% solve for the minimum 'a' (5 values) that will result in a minimum
rmse_ratio = @(a) rmse_id(a)/rmse_ref;
a = fmincon(rmse_ratio, ones([1,14]));
%Define final TF
TF = TF_temp(a);

end

function normVec = normVector(origMatrix)
    % calculate the distance vector of a given matrix, into a vector
	
    len = size(origMatrix, 1);
    normVec = zeros(len, 1);
    for x = 1:len
        normVec(x) = norm(origMatrix(x, :));
    end
end