%% SYSTEM SETTING

% In this code, we are trying to implement the concepts of level crossing
% rate (LCR) and average fade duration (AFD) of a randomly generated signal

%% CODE IMPLEMENTATION

clc;
clear all;
close all;

% Defining the necessary parameters
fD = 10;                        % maximum Doppler frequency in Hz
T = 10;                         % total time duration of the simulation
fs = 1000;                      % sampling frequency in Hz
N = fs*T;                       % number of samples
t = (0:N-1)/fs;                 % time vector

% Generating a Rayleigh fading signal (complex Gaussian Process)
rI = randn(1, N);
rQ = randn(1, N);
r_Rayleigh = sqrt(rI.^2 + rQ.^2); % Magnitude of the complex signal

% Define a threshold level Z for level crossing rate (LCR) and AFD
Z = 0.5;

%% Analytical Results

% Compute the level crossing rate (LCR)
dt = 1/fs;
derivative = diff(r_Rayleigh)/dt; % Derivative of the signal

% Find indices where the signal crosses the level Z downward
crossing_indices = find(r_Rayleigh(1:end-1) > Z & r_Rayleigh(2:end) <= Z & derivative(1:end) < 0);

% Level Crossing Rate (LCR)
LCR_analytical = length(crossing_indices) / T; % Crossings per second
disp(['Analytical Level Crossing Rate (LCR): ', num2str(LCR_analytical), ' crossings/second']);

% Calculate the Average Fade Duration (AFD)
% Find the time duration the signal stays below Z
below_threshold = r_Rayleigh < Z; % Binary vector where signal < Z
fade_durations = diff([0, below_threshold, 0]); % Find start and end of fades
fade_start_indices = find(fade_durations == 1); % Indices where fade starts
fade_end_indices = find(fade_durations == -1);  % Indices where fade ends

% Calculate the duration of each fade
fade_durations_sec = (fade_end_indices - fade_start_indices) / fs;

% Average Fade Duration (AFD)
AFD_analytical = mean(fade_durations_sec); % Average time spent below Z
disp(['Analytical Average Fade Duration (AFD): ', num2str(AFD_analytical), ' seconds']);

%% Using MATLAB Built-in or Utility Functions

% Use 'crossing' function to find level crossings
% crossing() returns the indices where the signal crosses a threshold

[zero_crossings, ~] = crossings(r_Rayleigh - Z);
LCR_builtin = length(zero_crossings) / T;
disp(['Built-in Level Crossing Rate (LCR): ', num2str(LCR_builtin), ' crossings/second']);

% Calculate Average Fade Duration using the crossing points
% Find segments where the signal stays below the threshold Z
fade_indices_builtin = r_Rayleigh < Z;
fade_diff_builtin = diff([0 fade_indices_builtin 0]);

fade_starts_builtin = find(fade_diff_builtin == 1);
fade_ends_builtin = find(fade_diff_builtin == -1);
fade_durations_builtin = (fade_ends_builtin - fade_starts_builtin) / fs;

% Built-in Average Fade Duration
AFD_builtin = mean(fade_durations_builtin);
disp(['Built-in Average Fade Duration (AFD): ', num2str(AFD_builtin), ' seconds']);

%% Compare Analytical and Built-in Results

disp('Comparison of Analytical and Built-in Results:');
disp(['Analytical LCR: ', num2str(LCR_analytical)]);
disp(['Built-in LCR: ', num2str(LCR_builtin)]);

disp(['Analytical AFD: ', num2str(AFD_analytical)]);
disp(['Built-in AFD: ', num2str(AFD_builtin)]);

% Plot the results
figure;
subplot(2,1,1);
plot(t, r_Rayleigh);
hold on;
yline(Z, '--r', 'Threshold Z');
xlabel('Time (s)');
ylabel('Signal Envelope');
title('Rayleigh Fading Signal with Threshold');

% Highlight the level crossings on the plot
plot(t(crossing_indices), r_Rayleigh(crossing_indices), 'ro', 'MarkerFaceColor', 'r');
legend('Signal', 'Threshold', 'Level Crossings');

subplot(2,1,2);
histogram(fade_durations_sec, 20);
xlabel('Fade Duration (s)');
ylabel('Count');
title('Histogram of Fade Durations');
