% Script that calcuates and plots MSE for differnt frame lengths

clear all; close all; clc;

% Read required file
[x_clean, Fs] = audioread(['Sound files/f1lcapae.wav']);
x_clean = x_clean(1:6.5e4);
noise = 0.01*randn(length(x_clean),1);
x_in = x_clean + noise;

% Compute estimate for PS of noise                    % Extract silent portion of data
Sn = var(noise).*length(x_in);             % Compute its power spectrum

% Pre aloocate y_out for speed
x_in = x_in'; y_out=0*x_in;

% Overlap parameters loop, testing frame length                                                    
j = 1;
for N = 128:4:1028;
    overlap = N/2;
    x=buffer(x_in,N,overlap);               % Place into matrix
    [N_samps,N_frames]=size(x);             % Number of samples and frames
    x_w=repmat(hanning(N),1,N_frames).*x;   % Window each frame
    
    clear X_w Y_w y_w_padded y_w;
    
    % Main loop for filtering
    for frame_no=1:N_frames-2
        
        % Compute FFTs
        X_w(:,frame_no)=fft([x_w(:,frame_no); zeros(1e3,1)]);             % Zero padd to make filtering more precise, avoidng digital distorsion
        Y_w(:,frame_no)=X_w(:,frame_no);
        
        % Compute filter gain, basic Weiner method
        snr = (abs(Y_w(:,frame_no)).^2-Sn)./Sn;                           % Signal to noise ratio
        for i = 1:length(snr)
            if snr(i) > 0
                filter_gain(i) = snr(i)/(1 + snr(i));
            else
                filter_gain(i) = 0;
            end;
        end
        
        % Filter signal
        Y_w = filter_gain'.*Y_w;
        
        % Signal back into time domain
        y_w_padded(:,frame_no)=ifft(Y_w(:,frame_no));
        y_w(:,frame_no) = y_w_padded(1:end-1e3,frame_no);                 % Remove zeros from zero padding
        
        % Add each filter frame to end of output
        y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=...
            y_out((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+...
            y_w(:,frame_no)';
    end
    
    % Calculate MSE
    MSE(j) = 1/length(x_in)*sum((x_clean'-y_out).^2);
    
    % Store MSE in array using index j
    if j == 1
        N_label = N;
    else
         N_label = [N_label, N];
    end
    
    j = j + 1;
end

% Plot results
figure; plot(N_label, MSE); 
xlabel('Frame Length'); ylabel('MSE');