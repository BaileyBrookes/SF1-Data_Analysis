clc; clear all; close all;

% Read data, exract induviual sounds
[data,Fs] = audioread('Sound files/piano_clean.wav');
n1 = data(1:3000);
n2 = data(3500:6000);
n3 = data(6001:8000);
n4 = data(8001:11760);
n5 = data(11761:14780);
n6 = data(14780:end);
notes = {data, n1, n2, n3, n4, n5, n6};
titles ={'Note 1', 'Note 2', 'Note 3', 'Note 4', 'Note 5', 'Note 6' }

% FFT of orginal data
Y = fft(data.*hann(length(data)));
L = length(Y);
Y = Y(1:L/2+1);
Y = abs(Y)/max(abs(Y));

% Normalised freq
f = Fs*(0:(L/2))/L;
[v,index] = min(abs(f-2e3));
f_norm = f*2*pi/Fs;

% Plot fft
figure;
plot(f(1:index-1),Y(1:index-1));
    
% Plot FFT for each note
figure; n = 1;
for i = 2:length(notes)
    % Window signal
    sig = notes{i}.*hann(length(notes{i}));
    
    % Compute FFT
    Y = fft(sig);
    L = length(Y);
    Y = Y(1:L/2+1);
    Y = abs(Y)/max(abs(Y));

    
    % Normalised freq
    f = Fs*(0:(L/2))/L;
    [v,index] = min(abs(f-2e3));
    f_norm = f*2*pi/Fs;
       
    % Plot fft
    subplot(3,2,n)
    plot(f(1:index-1),Y(1:index-1))
    xlabel('Frequency (Hz)')
    ylabel('Normalised Amplitude')
    title(titles{n})
    n = n+1;
end
