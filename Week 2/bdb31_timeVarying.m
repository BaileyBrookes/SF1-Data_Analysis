clear all; close all; clc;

% Read audio data
[x_in, Fs] = audioread('Sound files/dipper.wav');

% Pre aloocate y_out for speed
x_in = x_in'; y_out=0*x_in;

% Overlap parameters
N=512/2;                                % Frame size                                                                
overlap=256/2;                          % Overlap between frames                                                          
x=buffer(x_in,N,overlap);               % Place into matrix                                            
[N_samps,N_frames]=size(x);             % Number of samples and frames                             
x_w=repmat(hanning(N),1,N_frames).*x;   % Window each frame

% Compute estimate for PS of noise
noise = x(1:4e4);                     % Extract silent portion of data
Sn = var(noise).*length(x);             % Compute its power spectrum

% Main loop for filtering
for frame_no=1:N_frames-2
    
    % Compute FFTs
    X_w(:,frame_no)=fft([x_w(:,frame_no); zeros(1e3,1)]);             % Zero padd to make filtering more precise, avoidng digital distorsion
    Y_w(:,frame_no)=X_w(:,frame_no);
    
%     x_smooth = smoothdata(x_w(:,frame_no));
%     noise = x_w(:,frame_no) - x_smooth;
%     Sn = var(noise)*N;
    
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

% Plot results
subplot(2,1,1); plot(x_in); 
title('InputSignal'); xlabel('Sample Number');
subplot(2,1,2); plot(y_out);
title('Filtered Signal'); xlabel('Sample Number');

% Make player files with data in
original = audioplayer(x_in,Fs);
filtered = audioplayer(y_out,Fs);

% Ask to play results
disp('Press any key to here orgianal file'); waitforbuttonpress;
play(original);
disp('Please wait till audio stops playing and then press any key to here filtered file');
waitforbuttonpress;
play(filtered);

% Calculate and display MSE
MSE = 1/length(x_in)*sum((x_in-y_out).^2);
disp(['MSE is: ', num2str(MSE)])

% Wrtie the filtered files
audiowrite(filename,y_out,Fs)