clc; clear all; close all;

% Signal parameters
Fs = 100000;   % Sampling Frequency                            
T = 1/Fs;      % Period            
L = 150000;    % Length of signal         
t = (0:L-1)*T; % Time of samples
i = 1;         % Index for loop

for L = 1000:1000:10000
    
    % Generate Signal
    t = (0:L-1)*T;
    S = sin(2*pi*6000*t);
    W = hann(length(S));
    N = randn(1,length(S));
    
    % Window Signal
    X = W'.*S;
    X = [X,zeros(1,1e3)];
    
    % Compute one sided fft
    Y = fft(X);
    P2 = abs(Y)/max(abs(Y));
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    P1 = log(P1);

    % Compute Normalised Frequency 
    f = Fs*(0:(L/2))/L;
    f_norm = f*2*pi/Fs;

    % Plot fft
    hold on;
    plot(f_norm,P1)
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('Normalised frequency/rads^{-1}')
    ylabel('Normalised Amplitude')
    legendinfo{i} = [' N = ' num2str(L)]; 
    i = i + 1;
end
legend(legendinfo)