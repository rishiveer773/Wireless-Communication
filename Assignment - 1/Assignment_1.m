%% Wireless Communications Assignment-1 Report
%
%
%% SYSTEM SETTING
%
% In the below code, we will try to implement a system in which the input
% is an AWGN and we know the auto-correlation function of the received
% in-phase and quadrature components. Based on this information we have to
% now find the output function which is as follows:
%
% $$ r(t) = r_I(t)\cos(2\pi f_ct) - r_Q(t)\sin(2\pi f_ct) $$
%
% The auto-correlation of the channel modified input AWGNs is given as
% follows:
% $$ A_{r_I}(t, t+\tau) = \frac{P_r}{N_0} \sum_{i=0}^{N-1} \cos \left(
% \frac{2\pi v\tau}{\lambda} \cos(i\Delta\theta) \right) $$


%% IMPLEMENTATION

% TEST PUSH

clear all
close all
clc

Pr = 1;                                 % total input power received at the receiver
N0 = 10;                                % power of the AWGN input signal
X_f = N0 / 2;                           % Fourier transform of the AWGN
v = 50;                                 % velocity of the receiver
fc = 1e9;                               % carrier frequency of the input signal in Hz
c = 3 * 1e8;                            % speed of electromagnetic waves
angle = pi / 4;                         % angle at which the receiver is receiving the transmitted signal
tau = 0:0.001:1;                        % defining the time axis
f = linspace(-500, 500, length(tau));   % defining the frequency axis

% defining the input AWGN signals
x_I = randn(1, length(tau)) * sqrt(X_f);
x_Q = randn(1, length(tau)) * sqrt(X_f);

% computing the Doppler shift in the frequency caused due to the mobility
% of the receiver
fD = fc * v * cos(angle) / c ;

% we know the auto-correlation functions of the in-phase and quadrature
% components in terms of the Bessel function as follows
A_rI = Pr * besselj(0, 2*pi*tau*fD);    % auto-correlation function of the in-phase component
A_rQ = Pr * besselj(0, 2*pi*tau*fD);    % auto-correlation fucntion of the quadrature component

% hence, we can now find the PSD of the r_I(t) and r_Q(t) as it is
% equivalent to the Fourier transform of their respective auto-correlation
% functions
psd_rI = abs(fftshift(fft(A_rI)));
psd_rQ = abs(fftshift(fft(A_rQ)));

% thus, we can now find the equation of the frequency response of the
% respective in-phase and quadrature channels as follows
H_rI = sqrt(psd_rI) / X_f;
H_rQ = sqrt(psd_rQ) / X_f;

% finding the inverse Fourier transform of the obtained frequency response
h_rI = ifft(ifftshift(H_rI));
h_rQ = ifft(ifftshift(H_rQ));

% finding the actual r_I(t) and r_Q(t) now
rI = conv(x_I, h_rI);
rQ = conv(x_Q, h_rQ);
received_rI = rI(1:length(tau));
received_rQ = rQ(1:length(tau));

% finding the net output signal
r = received_rI .* cos(2*pi*fc*tau) - received_rQ .* sin(2*pi*fc*tau);

%% FIGURES

% plotting the input AWGN
figure
sgtitle('Input AWGN')
subplot(2,1,1)
plot(tau, x_I)
xlabel('Time')
ylabel('In-phase input AWGN')
grid on
subplot(2,1,2)
plot(tau, x_Q)
xlabel('Time')
ylabel('Quadrature input AWGN')
grid on

% plotting the auto-correlation functions of the received in-phase and
% quadrature components
figure
sgtitle('Auto-correlation functions of received signals')
subplot(2,1,1)
plot(tau, A_rI)
xlabel('Time')
ylabel('Auto-correlation function of in-phase component') 
grid on
subplot(2,1,2)
plot(tau, A_rQ)
xlabel('Time')
ylabel('Auto-correlation function of quadrature component')
grid on

% plotting the power spectral densities of the received in-phase and
% quadrature components
figure
sgtitle('PSDs of received signals')
subplot(2,1,1)
plot(f, psd_rI)
xlabel('Frequency')
ylabel('PSD of in-phase component')
grid on
subplot(2,1,2)
plot(f, psd_rQ)
xlabel('Frequency')
ylabel('PSD of quadrature component')
grid on

% plotting the frequency responses of the respective channels
figure
sgtitle('Frequency response of channel')
subplot(2,1,1)
plot(f, H_rI);
xlabel('Frequeny')
ylabel('Frequncy response of the in-phase channel')
grid on
subplot(2,1,2)
plot(f, H_rQ)
xlabel('Frequency')
ylabel('Frequency response of the quadrature channel')
grid on

% plotting the impulse responses of the respective channels
figure
sgtitle('Impulse response of channel')
subplot(2,1,1)
plot(tau, h_rI)
xlabel('Time')
ylabel('Impulse response of the in-phase channel')
grid on
subplot(2,1,2)
plot(tau, h_rQ)
xlabel('Time')
ylabel('Impulse response of the quadrature component')
grid on

% plotting the net received in-phase and quadrature signals
figure
sgtitle('Received in-phase and quadrature signals')
subplot(2,1,1)
plot(tau, received_rI)
xlabel('Time')
ylabel('Received in-phase signal')
grid on
subplot(2,1,2)
plot(tau, received_rQ)
xlabel('Time')
ylabel('Received quadrature component')
grid on

% plotting the net received signal
figure
sgtitle('Net received signal')
plot(tau, r)
xlabel('Time')
ylabel('Received signal')
grid on