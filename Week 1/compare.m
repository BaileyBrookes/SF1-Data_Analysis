clc; clear all; close all;

% Read audio file
[data,Fs] = audioread('Sound Files/armst_37_orig.wav');

% Main loop for computin dft and fft for differnt legnths N
time_index = 1;
for N = 1:10:length(data)
    data1 = data(1:N);

    % DFT
    tic;
    discrete = dft(data1,N);
    time_discrete(time_index) = toc;
    
    % FFT
    tic;
    fast = fft(data1);
    time_fast(time_index) = toc;
    
    time_index = time_index + 1;
end

% Plot results
figure; hold on;
plot(log(time_discrete))
plot(log(time_fast))
xlabel("Length of  transfrom, N")
ylabel("Log time to compute transfrom")
legend("DFT", "FFT")