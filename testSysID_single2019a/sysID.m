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

%Find the min and max 
min_val = min(sugar_vec);
max_val = max(sugar_vec);
%Find steady state value
steady = sugar_vec(end);

%In the example we ignore the peaks, min and max and just output a first order 
%response shifted by 160 
%Produce first order transfer function as output
s = tf('s');
% TF = -(max_val-min_val)*(4/(600*s+4)); %this one is like the reference
% ^^ just changed the value for steady state and simplified the transfer
% func

% % % Rise Time Tr
% % % : the time required for the waveform to go
% % % from 0.1 to 0.9 of the final value
% % % I Peak Time Tp: the time required to reach the first, or
% % % maximum, peak.
% % % I Peak value Mp: The magnitude value at the first, or
% % % maximum, peak.
% % % I Percent Overshoot, %OS: The amount that the waveform
% % % overshoots the steady-state, or final value at the peak
% % % value, expressed as a percentage of the steady-state value.
% % % I Settling Time Ts
% % % : the time required for the transient’s
% % % damped oscillations to reach and stay within ±2% of the
% % % steady-state value


% OS
OS = abs((minPKS(1)-steady)/steady);

%time of first peak (in our case a local minima)
tp = minLOCS(1);

%ts = +-2%
ts = 0;
tsarray = logical((steady-0.02*steady)<sugar_vec & sugar_vec<(steady+0.02*steady));
for i = 0:length(tsarray)-1
    i_rev = length(tsarray)-i;
    if(tsarray(i_rev) == 0)
        if (i_rev>ts)
            ts = i_rev;
        end
    end
end

val90 = (steady) +(max_val - steady)*0.1;
val10 = (steady) +(max_val - steady)*0.9;
array90 = logical(sugar_vec<val90);
array10 = logical(sugar_vec<val10);
t10 = find(array10, 1, 'first');
t90 = find(array90, 1, 'first');
tr = (t90-t10)*1.6;

damp = abs(log(OS)/sqrt( pi^2 + log(OS)^2 )) * 0.84;
% wn = (pi/tp)/sqrt(1-damp); 
%failed attempt
wn = (pi/(tp*1.1)) / (sqrt(1-damp^2));
%when the first minimum is less than the steady state (i.e. it did the
%opposite of overshoot)

TF = (wn^2 + 10^(-10)) / (( (s^2) + (2*wn*damp*s) + (wn^2))); 

TF = TF * -(max_val - steady);

if (minPKS(1)-1) > steady
    %new pole
     p = 1.75;
    damp = 0.95;
    wn = ( (pi/(tp)) / (sqrt(1-damp^2)) )*0.43;
    
   
    TF = (wn^2) / (( (s^2) + (2*wn*damp*s) + (wn^2))*(s+p)); 
    TF = TF * -1.745*(max_val - steady);
    a = 0
    
elseif length(maxPKS)>1
    
    if maxPKS(2)>(steady+1)
%         new pole
        p = 1.;
        damp = 0.7;
        wn = ((pi/tp) / (sqrt(1-damp^2)))*0.9;
        
        TF = (wn^2) / (( (s^2) + (2*wn*damp*s) + (wn^2))*(s+p)); 
        TF = TF * -1*(max_val - steady);
        a = 1
    end
elseif length(minPKS) < 2
        damp = 1;
        wn = ((pi/tp) / (sqrt(1-damp^2)))*1.2;
        
        TF = (wn^2) / (( (s^2) + (2*wn*damp*s) + (wn^2))); 
        TF = TF * -1*(max_val - steady);
        a = 2

end



% wn = 1.8/tr;
% damp = 4/(wn*ts)

% % stepinfo(sugar_vec, time_vec);

% wn = 1.8/(minLOCS(1)/2.5) %~0.0054
% damp = 0.9;
% % %^^ works okay for some things (better than ref) 
% % % except when second oscillation is quite large
% % %or if steady state difference is small


%Produce initial condition (offset from zero)
IC = max_val; %changed the start value to the max value

end