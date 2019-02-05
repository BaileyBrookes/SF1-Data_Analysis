% Sript that filters 1 channel audio files uisng the basic Weiner filter
% method, predicitng noise  by subtracting the smoothed versionof a frame
% from its orginal


clear all; close all; clc;

% Ask user for audio file to be filtered?
disp('Which audio file would you like to filter? Type the exact file name')
disp('with the .wav extenstion, ensure its in the same folder as this script');
audio_file = input('and inclued the quote marks: ');

% Read required file
[x_in, Fs] = audioread(audio_file);
x_in = x_in(:,1)'; y_out=0*x_in;

% Overlap parameters
N=256;                                % Frame size                                                                
overlap=128;                          % Overlap between frames                                                          
x=buffer(x_in,N,overlap);               % Place into matrix                                            
[N_samps,N_frames]=size(x);             % Number of samples and frames                             
x_w=repmat(hanning(N),1,N_frames).*x;   % Window each frame

% Main loop for filtering
for frame_no=1:N_frames-2
    
    % Compute FFTs
    X_w(:,frame_no)=fft([x_w(:,frame_no); zeros(1e3,1)]);             % Zero padd to make filtering more precise, avoidng digital distorsion
    Y_w(:,frame_no)=X_w(:,frame_no);
    
    % Another attempt at noise
    x = x_w(:,frame_no);
    x_smooth = smoothdata(x);
    noise = x - x_smooth;
    Sn = var(noise)*length(x);
    
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
subplot(1,1,1); plot(x_in); 
title('Input Signal'); xlabel('Sample Number');
subplot(2,1,2); plot(y_out);
title('Filtered Signal'); xlabel('Sample Number');

% Write the filtered files
filename = [audio_file, '_filtered_smooth.wav'];
audiowrite(filename,y_out,Fs)

% Calculate and display MSE
MSE = 1/length(x_in)*sum((x_in-y_out).^2);
disp(['MSE is: ', num2str(MSE)])