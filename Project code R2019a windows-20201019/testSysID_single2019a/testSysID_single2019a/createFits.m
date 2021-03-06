function [fitresult, gof] = createFits(time_vec, patient_sugar_resp)
%CREATEFITS(TIME_VEC,PATIENT_SUGAR_RESP)
%  Create fits.
%
%  Data for 'untitled fit 1' fit:
%      X Input : time_vec
%      Y Output: patient_sugar_resp
%  Data for 'untitled fit 2' fit:
%      X Input : time_vec
%      Y Output: patient_sugar_resp
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Oct-2020 20:00:58

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 2, 1 );
gof = struct( 'sse', cell( 2, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time_vec, patient_sugar_resp );

% Set up fittype and options.
ft = fittype( 'sin2' );
excludedPoints = excludedata( xData, yData, 'Indices', 1 );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf];
opts.Robust = 'LAR';
opts.StartPoint = [172.975388624955 0.00218317766059054 0.0701418736331849 77.9421263722672 0.00436635532118109 1.65656668841431];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult{1}, xData, yData, excludedPoints );
legend( h, 'patient_sugar_resp vs. time_vec', 'Excluded patient_sugar_resp vs. time_vec', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'time_vec', 'Interpreter', 'none' );
ylabel( 'patient_sugar_resp', 'Interpreter', 'none' );
grid on

%% Fit: 'untitled fit 2'.
[xData, yData] = prepareCurveData( time_vec, patient_sugar_resp );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );

% Fit model to data.
[fitresult{2}, gof(2)] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'untitled fit 2' );
h = plot( fitresult{2}, xData, yData );
legend( h, 'patient_sugar_resp vs. time_vec', 'untitled fit 2', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'time_vec', 'Interpreter', 'none' );
ylabel( 'patient_sugar_resp', 'Interpreter', 'none' );
grid on


