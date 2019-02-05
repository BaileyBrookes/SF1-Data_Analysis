clear all; close all; clc;

f = 100;
Fs = 1000;
T = 1/Fs;
t = 0:T:1;
L = length(t);

sig = sin(2*pi*f*t);

spec = log(abs(fft(sig)));
spec = spec(1:round(L/2+1));
f = Fs*(0:(L/2))/L;

plot(f,spec(1:length(f)))