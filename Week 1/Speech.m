clc; clear all; close all;
% Read data, extract induviual sounds
[data,Fs] = audioread('Sound files/f1lcapae.wav');
data = data(1:3E4);
d1 = data(1:5390);
d2 = data(5391:1.117e4);
d3 = data(1.1171e4:1.483e4); 
d4 = data(1.4831e4:1.639e4);
d5 = data(1.639e4:1.870e4);
d6 = data(1.870e4:2.097e4);
d7 = data(2.097e4:2.3e4);
d8 = data(2.3e4:2.485e4);
d9 = data(2.485e4:end);
sounds = {data, d1, d2, d3, d4, d5, d6, d7, d8, d9};
ma = {'Entire phrase', 'Sound 1', 'Sound 2', 'Sound 3', 'Sound 4', 'Sound 5', 'Sound 6', 'Sound 7', 'Sound 8', 'Sound 9'};

% Window each one
win = 'hann'; 
for i = 2:length(sounds)
    w = window(str2func(win), length(sounds{i}));
    sounds{i} = sounds{i}.*w;
end

% Loop for ploting results
n = 1;
for i = 1:length(sounds)
    % Compute FFT
    Y = fft([sounds{i};zeros(1e6,1)]);
    L = length(Y);
    Y = Y(1:round(L/2)+1);
    Y = abs(Y)/max(abs(Y));
    
    % Normalised freq
    f = Fs*(0:(L/2))/L;
    [v,index] = min(abs(f-1e3));
    f_norm = f*2*pi/Fs;
    
    Yfft{i} = [f',Y(1:length(f))];
    
    % Plot fft
    subplot(5,2,n)
    plot(f(1:index-1),Y(1:index-1))
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title(ma{n})
    n = n+1;
    sound(sounds{i})
    waitforbuttonpress
end
 




