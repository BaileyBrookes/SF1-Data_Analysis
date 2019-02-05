clc; clear all;

% Specify parameters
Fs = 10e3;
T = 1/Fs;
L = 1000;
t = (0:L-1)*T;
i = 1;
windows = {'rectwin','hann', 'hamming', 'triang'};

% Generate Signal
t = (0:L-1)*T;
S = sin(2*pi*60*t);
N = randn(1,length(S));

% Normalised freq
f = Fs*(0:(L/2))/L;
f_norm = f*2*pi/Fs;

% Window signal
W = window(str2func(windows{4}),length(S));
X = W'.*S; %+N;
X = [X,zeros(1,1e3)];

% Compute one sided fft
Y = fft(X);
P1 = Y(1:L/2+1);
P1 = abs(P1)/max(abs(P1));
P1 = log(P1);
data{i} = P1;

% Plot fft
hold on;
plot(f_norm,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Normalised frequency/rads^{-1}')
ylabel('|FFT|')