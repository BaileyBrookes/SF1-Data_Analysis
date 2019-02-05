clc; clear all; close all;
% Read data
[data,Fs] = audioread('Sound files/organ.wav');
T = 1/Fs;

% Plot data in time domain 
plot(data(:,1))
xlabel('Sample')
ylabel('Amplitude')

% Window and zero padd data
X = [data(:,1).*hann(length(data)); zeros(1e6,0)];

% Compute FFT of data
Y = fft(X);
L = length(Y);
Y = Y(1:L/2+1);
Y = abs(Y)/max(abs(Y));

% Normalised freq
f = Fs*(0:(L/2))/L;
[v,index] = min(abs(f-3e3));
f_norm = f*2*pi/Fs;

% Plot fft
figure;
plot(f(1:index-1),Y(1:index-1));
xlabel('Frequency (Hz)')
ylabel('Noralised Amplitude')
