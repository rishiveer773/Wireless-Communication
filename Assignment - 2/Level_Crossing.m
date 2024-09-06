%% SYSTEM SETTING

% In this code, we are trying to implement the concepts of level crossing
% rate (LCR) and average fade duration (AFD) of a randomly generated signal


%% CODE IMPLEMENTATION

clc;
clear all;
close all;

% Defining the necessary parameters
fD = 240;                                       % maximum Doppler frequency in Hz
T = 10;                                         % total time duration of the simulation
fs = 1000;                                      % sampling frequency in Hz
N = fs*T;                                       % number of samples
t = (0:N-1)/fs;                                 % time vector

% Generating a Rayleigh fading signal (complex Gaussian Process)
rI = randn(1, N);
rQ = randn(1, N);
r_Rayleigh = sqrt(rI.^2 + rQ.^2);

% Define a threshold level Z for level crossing rate (LCR) and AFD
Z = 1;

% Compute the level crossing rate (LCR)
dt = 1/fs;
derivative = diff(r_Rayleigh)/dt;               % derivative of the signal

% Finding the indices where the signal crosses the level Z downward
crossing_indices = find(r_Rayleigh(1:end-1) > Z & r_Rayleigh(2:end) <= Z & derivative(1:end) < 0);
LCR = length(crossing_indices) / T;             % Level Crosssing Rate
disp(['Level Crossing Rate (LCR): ', num2str(LCR), ' crossings/second']);

% Calculating the average fade duration. First, let us find the time 
% duration for whihc the signal stays below Z
below_threshold = r_Rayleigh < Z;               % Binary vector where signal < Z
fade_durations = diff([0, below_threshold, 0]); % Find start and end of fades
fade_start_indices = find(fade_durations == 1); % Indices where fade starts
fade_end_indices = find(fade_durations == -1);  % Indices where fade ends

% Calculating the duration of each fade
fade_durations_sec = (fade_end_indices - fade_start_indices) / fs;

% Average Fade Duration (AFD)
AFD = mean(fade_durations_sec);                 % average time spent below Z
disp(['Average Fade Duration (AFD): ', num2str(AFD), ' seconds']);

% Checking with the formula
rho = Z / sqrt(2);                              % sqrt(2) because variance is 1
formula_lcr = sqrt(2 * pi) * fD * rho * exp(-(rho^2));
disp(['Level Crossing Rate (LCR) by Formula: ', num2str(formula_lcr), ' crossings/second'])

avg_fade_duration = (exp(rho^2) - 1) / (rho * fD * sqrt(2*pi));
disp(['Average Fade Duration (AFD) by Formula: ', num2str(avg_fade_duration), ' seconds']);


%% FIGURES

figure;
subplot(2,1,1);
plot(t, r_Rayleigh);
hold on;
yline(Z, '--r', 'Threshold Z');
xlabel('Time (s)');
ylabel('Signal Envelope');
title('Rayleigh Fading Signal with Threshold');

% Highlight the level crossings on the plot
% plot(t(crossing_indices), r_Rayleigh(crossing_indices), 'ro', 'MarkerFaceColor', 'r');
% legend('Signal', 'Threshold', 'Level Crossings');

subplot(2,1,2);
histogram(fade_durations_sec, 20);
xlabel('Fade Duration (s)');
ylabel('Count');
title('Histogram of Fade Durations');
