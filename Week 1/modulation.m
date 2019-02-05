clc; clear all; close all;

%% Specify parameters
Fs = 100e3; 
T = 1/Fs;
L = 1000;

% Generate Signal
t = (0:L-1)*T;
S = sin(2*pi*4000*t);

% Normalised freq
f = Fs*(0:(L/2))/L;
f_norm = f*2*pi/Fs;

%% Random Modulation scheme, has same effect as noise
len = length(S);
a1 = zeros(1,len);
a1(1) = 0;
for i = 2:len
    a1(i) =a1(i-1) + 0.1*randn(1); 
end
X1 = a1.*S;
% X1 = [X1,zeros(1,1e3)];

%% Linear Increase, smoothing effect
A = 1; B = 1;
a2(1) = 0;
for i = 2:len
    a2(i) =A*a2(i-1) + B; 
end
X2 = a2.*S;
% X2 = [X2,zeros(1,1e3)];

%% Periodic modulation, like standard AM
beta = 0.5; phase = 0.01;
for i = 1:len
    a3(i) = 1 + beta*sin(2*pi*phase*i); 
end
X3 = a3.*S;   
% X3 = [X3,zeros(1,1e3)]; 

scheme = {S, X1, X2, X3};

%% Compute one sided fft, choosing a modulation scheme
for i = 1:length(scheme)
    Y = fft(scheme{i});
    P1 = Y(1:L/2+1);
    P1 = abs(P1)/max(abs(P1));
    P1 = log(P1);
    data{i} = P1;
end

%% Plot fft
figure; hold on;
for i = 1:length(scheme)
    plot(f_norm, data{i})
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('Normalised frequency/rads^{-1}')
    ylabel('Normalised Amplitude')
end