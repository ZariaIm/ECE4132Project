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

%%  system identification
% %% polynomial fit
% 
% % set useful params
% IC = sugar_vec(1);
% % extend data with constant final value, so polyfit smoothens at the end
% time = 1:2000;
% order = 5;
% sugar_vec_extended = [sugar_vec, sugar_vec(end)*ones(1,time(end)-time_vec(end)-1)];
% coeffs = polyfit(time, sugar_vec_extended-IC, order);   % fit curve of order n to patient data
% 
% % add step response for each term of the polynomial to TF
% TF = 0;
% s = tf('s');
% for i = 1:length(coeffs)
%     s_order = order - i + 1;
%     TF = TF + coeffs(i)*factorial(s_order)/s^s_order; 
% end
% TF = minreal(TF)
% 
% [z, p, k] = zpkdata(TF);
% z = cell2mat(z)
% p = cell2mat(p)

% optimisation
% Produce initial condition (offset from zero)
IC = sugar_vec(1); %changed the start value to the max value
% simulate the reference as comparison
[TF_ref, IC_ref] = referenceID(patient);
Y_ref = step(TF_ref,time_vec);
ref_resp = Y_ref+IC_ref;
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
[a, ratio] = fmincon(rmse_ratio, ones([1,14]));
ratio;
%test what those minimum values are giving.
TF = TF_temp(a);
% [z, p, k] = zpkdata(TF);
% z = cell2mat(z)
% p = cell2mat(p)
rmse_id_test = rmseFct(sugar, id(a));


% %% old code
% if ratio>0.9
%         %Find the peaks in the data, this is when the slope changes sign
%     [maxPKS,maxLOCS] = findpeaks(sugar_vec,time_vec);
%     [minPKS,minLOCS] = findpeaks(-sugar_vec,time_vec);
%     minPKS = -minPKS; 
% 
%     %Find steady state value
%     steady = sugar_vec(end);
%     s = tf('s');
% 
%     if length(minPKS)<1
%         damp = 7.5;
%         wn = 2;
%         a = 4/570;
%     else
% 
%         % OS
%         OS = (steady-minPKS(1))/steady;
%         %time of first peak (in our case a local minima)
%         tp = minLOCS(1);
%             damp = abs(log(OS)/sqrt( pi^2 + log(OS)^2 )) * 0.8;
%             wn = (pi/tp) / (sqrt(1-damp^2))*0.85;    
%             a = wn*damp*steady*0.33;
% 
%         if OS < 0 || minLOCS(1)>600
%             damp = 7.5;
%             wn = 2;
%             a = 4/570;
%         elseif length(maxPKS)>1
%             if maxPKS(2)>(steady+1)
%             damp = 0.85;
%             wn = 0.75*(pi/tp) / (sqrt(1-damp^2));
%             a = wn*damp*steady;
%             end
%         end
%     end
%         %Transfer function that is used for all 'types' of systems
%         TF = -(sugar_vec(1) - steady)*a*(wn^2 + 10^(-10)) / ((s+a)*( (s^2) + (2*wn*damp*s) + (wn^2)));
% %calculate the rmse of the reference
% rmseFct = @(x, y) sqrt(sum((normVector(x - y)).^2)/(size(x, 1)));
% %make function handles for the id
% id =step(TF ,time_vec) + IC;
% 
% rmse_id = rmseFct(sugar, id);
% 
% %solve for the minimum 'a' (5 values) that will result in a minimum
% rev_rmse_ratio = rmse_id/rmse_ref
% 
% end
end