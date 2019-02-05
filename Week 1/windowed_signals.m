clc; clear all; close all;

Fs = 100000;                             
T = 1/Fs;                
L = 150000;             
t = (0:L-1)*T;  
i = 1;

for L = 1:1000:10000
    
    % Generate Signal
    t = (0:L-1)*T;
    S = sin(2*pi*6000*t);
    W = hann(length(S));
    N = randn(1,length(S));
    X = W'.*S;%+1*N;

    % figure;
    % plot(1000*t,X)
    % title('Signal Corrupted with Zero-Mean Random Noise')
    % xlabel('t (milliseconds)')
    % ylabel('X(t)')
    
    % Compute one sided fft
    Y = fft(X);
    P2 = abs(Y)/max(abs(Y));
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1 = log(P1);

    f = Fs*(0:(L/2))/L;
    f_norm = f*2*pi/Fs;

    % Plot fft
    hold on;
    plot(f_norm,P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('Normalised frequency/rads^{-1}')
    ylabel('|P1(f)|')
    legendinfo{i} = [' N = ' num2str(L)]; 
    i = i + 1;
end
legend(legendinfo)