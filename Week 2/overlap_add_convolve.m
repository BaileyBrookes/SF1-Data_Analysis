tic;
x_in=0.5*cos([1:10000]*pi/4)+sin([1:10000]*pi/100)+randn(1,10000);    % Generate Noisey signal
y_out_cov=0*x_in;                                                         % Pre allocate y for speed

N=512;                                                                % N = Number of datapoints
overlap=256;                                                          % Size of overlap = N/2

x=buffer(x_in,N,overlap);                                             % Split x into matrix, each with N rows and overlapping with size overlap

[N_samps,N_frames]=size(x);                                           % Get samples per frame (#rows) and number of frams (#cols)

x_w=repmat(hanning(N),1,N_frames).*x;                                 % Window samples with hann window, (aviods for loops)

% Main loop for filtering
for frame_no=1:N_frames-2    
    
    y_w(:,frame_no)=y_w(:,frame_no);                                  % Set output FFT to inputs FFT (pre allocate for speed)
    
    % Simple attempt at noise reduction, multipling parts of fft by fixed
    % gains
    y_w(2:N/8,frame_no)=conv(0.1,x_w(2:N/8,frame_no));                      % Set output 2 through N/8 to 1/10 of the input 
    y_w(N/4+1:N/2,frame_no)=conv(0.2,x_w(N/4+1:N/2,frame_no));              % Set output N/4+1 down to N/2 to 1/5 of input
    
   % Y_w(N:-1:N/2+2,frame_no)=conj(Y_w(2:N/2,frame_no));               % Set second half of frame to complex conjuage of first half. AS overlap is half this bit is filtered in the next stage
    
   % y_w(:,frame_no)=ifft(Y_w(:,frame_no));                            % Take iFFT
    
    y_out_cov((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)=...          % Construct output by adding each frame back togeather
        y_out_cov((frame_no-1)*overlap+1:(frame_no-1)*overlap+N)+...
        y_w(:,frame_no)';
    
end
time_conv = toc