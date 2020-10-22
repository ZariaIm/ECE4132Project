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
% % %Find the peaks in the data, this is when the slope changes sign
% % [maxPKS,maxLOCS] = findpeaks(sugar_vec,time_vec);
% % [minPKS,minLOCS] = findpeaks(-sugar_vec,time_vec);
% % minPKS = -minPKS; 
% % 
% % %Find the max 
% % max_val = max(sugar_vec);
% % %Find steady state value
% % steady = sugar_vec(end);
% % s = tf('s');
% % 
% % if length(minPKS)<1
% %     damp = 7.5;
% %     wn = 2;
% %     a = 4/570;
% % else
% %     
% %     % OS
% %     OS = (steady-minPKS(1))/steady;
% %     %time of first peak (in our case a local minima)
% %     tp = minLOCS(1);
% %         damp = abs(log(OS)/sqrt( pi^2 + log(OS)^2 )) * 0.8;
% %         wn = (pi/tp) / (sqrt(1-damp^2))*0.85;    
% %         a = wn*damp*steady*0.33;
% % 
% %     if OS < 0 || minLOCS(1)>600
% %         damp = 7.5;
% %         wn = 2;
% %         a = 4/570;
% %     elseif length(maxPKS)>1
% %         if maxPKS(2)>(steady+1)
% %         damp = 0.85;
% %         wn = 0.75*(pi/tp) / (sqrt(1-damp^2));
% %         a = wn*damp*steady;
% %         end
% %     end
% % end
% %     %Transfer function that is used for all 'types' of systems
% %     TF = -(max_val - steady)*a*(wn^2 + 10^(-10)) / ((s+a)*( (s^2) + (2*wn*damp*s) + (wn^2)));
% %     %Produce initial condition (offset from zero)
% %     IC = max_val; %changed the start value to the max value



% %    Coefficients (with 95% confidence bounds):
% %     Sum of 3 sines
% %        a1 =       145.3  (114, 176.6)
% %        b1 =    0.001319  (0.0006594, 0.001978)
% %        c1 =       1.056  (0.4372, 1.674)
% %        a2 =       167.8  (-1527, 1863)
% %        b2 =    0.004827  (0.001147, 0.008508)
% %        c2 =       1.546  (-1.558, 4.65)
% %        a3 =         130  (-1582, 1842)
% %        b3 =    0.005414  (0.00203, 0.008799)
% %        c3 =       4.188  (1.312, 7.063)
% %     Sum of 2 sines
% %        a1 =        3541  (-3.445e+06, 3.452e+06)
% %        b1 =    0.001067  (-0.04375, 0.04589)
% %        c1 =       2.013  (-30.91, 34.94)
% %        a2 =        3268  (-3.446e+06, 3.452e+06)
% %        b2 =    0.001156  (-0.04375, 0.04606)
% %        c2 =       5.091  (-26.71, 36.89)

% Range for IC should be from -200 to 200 at most
    s = tf('s');
%Use the following parameters and optimise for the lowest RMSE possible
 % parameters
   a1 = 3541;
   b1 = 0.001063;
   c1 = 2.012;
   a2 = 3262;
   b2 = 0.001156;
   c2 = 5.091; 
%    a3 = 130;
%    b3 = 0.005414;
%    c3 = 4.188;
   IC = 0;
   
   n1_1 = a1*sin(c1);
   n1_2 = a2*sin(c2);
%    n1_3 = a3*sin(c3);
   n2_1 = a1*b1*cos(c1);
   n2_2 = a2*b2*cos(c2);
%    n2_3 = a3*b3*cos(c3);
   d_1 = b1^2;
   d_2 = b2^2;
%    d_3 = b3^2;
   
%3 sines
% TF = ((n1_1*s^2 +n2_1*s)/(s^2 + d_1) + (n1_2*s^2 +n2_2*s)/(s^2 + d_2) + (n1_3*s^2 +n2_3*s)/(s^2 + d_3))

%2 sines
TF = ((n1_1*s^2 +n2_1*s)/(s^2 + d_1) + (n1_2*s^2 +n2_2*s)/(s^2 + d_2))
%MANUAL
TF = tf([170.9 -0.2135 0.0008539 -5.73e-07 0], [1 1e-06 2.43e-06 0 1.4e-12])
[p,z] = pzmap(TF)

% %     %find polynomial coefficients in time domain
% %     pfit = polyfit(sugar_vec, time_vec, 7);
% %     poly = pfit(1)*s^2 + pfit(2)*s^1 + pfit(3);
% %    ilaplace(poly)
% %     scale = (max(sugar_vec) - sugar_vec(end));
% %     TF = 1/(pfit(1)*s^2 + pfit(2)*s^1 + pfit(3));
% % syms x
% % % polynomial fit for patient response
% % f = 4.569*x^7 + -4.42*x^6 + -24.3*x^5 + 25.75*x^4 + 27.19*x^3 + -25.31*x^2 + -9.27*x + 95.07;
% % temp = laplace(f)
% % TF = 1.8*tf([-0.2545 0.0006536*2 -3.043*6*10^(-7) -1.364*10^(-9)*24 (2.415*10^(-12))*120 (-1.523*10^(-15))*720 3.41*10^(-19)*5040], [1 0 0 0 0 0 0 0])
% 

% pp = spline(time_vec,sugar_vec,time_vec);


% % % [xData, yData] = prepareCurveData( time_vec, sugar_vec );
% % % 
% % % % Set up fittype and options.
% % % ft = fittype( 'sin3' );
% % % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % % opts.Display = 'Off';
% % % opts.StartPoint = [147.440007136686 0.00218166156499291 0.266362399817917 76.6301159070064 0.00436332312998582 1.89008192518712 20.9143646972681 0.00872664625997165 1.7102165411476];
% % % 
% % % % Fit model to data.
% % % [fitresult, gof] = fit( xData, yData, ft, opts );
% % % 
% % % % Plot fit with data.
% % % figure( 'Name', 'untitled fit 3' );
% % % h = plot( fitresult, xData, yData );
% % % legend( h, 'sugar_vec vs. time_vec', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % % % Label axes
% % % xlabel( 'time_vec', 'Interpreter', 'none' );
% % % ylabel( 'sugar_vec', 'Interpreter', 'none' );
% % % grid on
% % % 
% % % 
% % % syms x
% % % s6c = coeffvalues(fitresult);
% % % a = 1:3:length(s6c);
% % % s6ca = s6c(a);
% % % b = 2:3:length(s6c);
% % % s6cb = s6c(b);
% % % c = 3:3:length(s6c);
% % % s6cc = s6c(c);
% % % for len = 1:(length(s6c)/3)
% % %     f(len) = s6ca(len)*sin(s6cb(len)*x+s6cc(len));
% % %     flap(len) = laplace(f(len));
% % % end
% % % TFtest = sum(flap);
% % % [num,den] = numden(TFtest);
% % % n = double(coeffs(vpa(num)))
% % % d = double(coeffs(vpa(den)))
% % % 
% % % 
% % % TF = tf(d, n);
% % % IC = 0;

% [xData, yData] = prepareCurveData( time_vec, sugar_vec );
% 
% % Set up fittype and options.
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.SmoothingParam = 7.70916786753862e-06;
% 
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % Plot fit with data.
% figure( 'Name', 'untitled fit 3' );
% h = plot( fitresult, xData, yData );
% legend( h, 'sugar_vec vs. time_vec', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'time_vec', 'Interpreter', 'none' );
% ylabel( 'sugar_vec', 'Interpreter', 'none' );
% grid on



% % [xData, yData] = prepareCurveData( time_vec, sugar_vec );
% % 
% % % Set up fittype and options.
% % ft = fittype( 'gauss3' );
% % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % opts.Display = 'Off';
% % opts.Lower = [-50 -50 0 -50 -50 0 -50 -50 0];
% % opts.StartPoint = [152.944541290625 0 158.651806870014 87.3714291413174 286 388.575052852056 80.963493807267 1440 409.708654307612];
% % 
% % % Fit model to data.
% % [fitresult, gof] = fit( xData, yData, ft, opts );
% % rmse = gof.rmse
% % % % Plot fit with data.
% % % figure( 'Name', 'untitled fit 3' );
% % % h = plot( fitresult, xData, yData );
% % % legend( h, 'sugar_vec vs. time_vec', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % % % Label axes
% % % xlabel( 'time_vec', 'Interpreter', 'none' );
% % % ylabel( 'sugar_vec', 'Interpreter', 'none' );
% % % grid on
% % c = coeffvalues(fitresult);
% % syms x
% % sys = c(1)*exp(-((x-c(2))/c(3))^2) + c(4)*exp(-((x-c(5))/c(6))^2) + c(7)*exp(-((x-c(8))/c(9))^2)
% % L = laplace(sys)
% % [n,d] = numden(L);
% % n = vpa(n,2)
% % n = coeffs(n)
% % n = double(n)
% % d = vpa(d,2)
% % d = coeffs(d)
% % d = double(d)
% % bode(tf(n,d))
% % 
% % 
% % % Set up fittype and options.
% % [xData, yData] = prepareCurveData( time_vec, sugar_vec );
% % ft = fittype( 'fourier3' );
% % opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % opts.Display = 'Off';
% % opts.StartPoint = [0 0 0 0 0 0 0 0.00218166156499291];
% % % % 
% % % Fit model to data.
% % [fitresult, gof] = fit( xData, yData, ft, opts );
% % gof.rmse
% % % % % Plot fit with data.
% % % % figure( 'Name', 'untitled fit 3' );
% % % % h = plot( fitresult, xData, yData );
% % % % legend( h, 'sugar_vec vs. time_vec', 'untitled fit 3', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % % % % Label axes
% % % % xlabel( 'time_vec', 'Interpreter', 'none' );
% % % % ylabel( 'sugar_vec', 'Interpreter', 'none' );
% % % % grid on
% % 
% % %time response
% % fc = coeffvalues(fitresult);
% % a0 = fc(1);
% % a1 = fc(2);
% % b1 = fc(3);
% % a2 = fc(4);
% % b2 = fc(5);
% % a3 = fc(6);
% % b3 = fc(7);
% % w = fc(8);
% % f = a0 + a1*cos(time_vec*w) + b1*sin(time_vec*w) + a2*cos(time_vec*w) + b2*sin(time_vec*w)+ a3*cos(time_vec*w) + b3*sin(time_vec*w);
% % 
% % % % TF = tf();
% % % % IC = 0;
% % 
% % 
% % %freq response
% % fs=12000; %sampling frequency (arbitary prolly there is something built in as well 
% % resp = fft(f); %fourier transform 
% % len = length(resp); 
% % 
% % 
% % s_number = pow2(nextpow2(len)); %returns sample points multiple of 2 %used for no apparent reason just to see 
% % y_s_number = fft(f,s_number); %fourier transform with desired sample points
% % freq = (0:s_number-1)*(fs/s_number)/10; % frequency vector
% % pow = abs(y_s_number).^2/s_number;   % power spectrum 
% % 
% % resp_shift = fftshift(resp); %fourier transform with shifted response
% % pow_shift = abs(resp_shift).^2/len; %power spec of shifter fourier transform
% % freq_shift = (-len/2:len/2-1)*(50/len); %freuency rnage or spectrum (50 is a random value)
% %    
% % 
% % %plot fourier with desired smaple points
% % figure(2)
% % plot(freq(1:floor(s_number/2)),pow(1:floor(s_number/2)))
% % xlabel('Freq')
% % ylabel('Pow')
% % title("fourier_pow2")
% % 
% % %plot fourier shift
% % figure(3)
% % plot(freq_shift,pow_shift)
% % xlabel('Freq')
% % ylabel('Pow')
% % title("Shifted fft")



end