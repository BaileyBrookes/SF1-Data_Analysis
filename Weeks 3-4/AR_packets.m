clear all; close all; clc;

% Read data
[x_in, Fs] = audioread('Sound Files/armst_37_orig.wav.'); 


% Set sections to zero
x_OG = x_in;
x_in(6500:6500+100) = 0;
x_in(8888:8888+100) = 0;
x_in(1.356E4:1.356E4+100) = 0;
data_length = length(x_in);
[index] = find(~x_in);

% Sepcify Model Order
P = 100;

% Error
m_e = 0; var_e=0;
for zero_location = [6500 8888 1.356E4];

    x = x_in(zero_location-500:zero_location+P+500);
    
    % Form generator matrix
    N = length(x); clear G;
    for i = 1:P
        G(:,i) = [zeros(i,1); x(1:N-i)];
    end
    
    % Form parpameter estimate fro MAP and ML
    m_theta = zeros(P,1); cov_theta = 5*eye(P);
    phi = G'*G + var_e*inv(cov_theta);
    The = G'*x  + var_e*inv(cov_theta)*m_theta;
    theta_MAP = inv(phi)*The;
    theta_MLP = inv(G'*G)*G'*x;
    
    x_MAP = G * theta_MAP;
    
    % No predict the empty bit
    a = linspace(1,0,100);
    x_seg_f = x_in(zero_location - P :zero_location-1);
    x_seg_b = x_in(zero_location + 101 :zero_location + 100 + P);
    
    for i = 1:100
        x_f(i) = theta_MAP'*x_seg_f(end-P+1:end);
        x_b(i) = theta_MAP'*x_seg_b(1:P);
        x_seg_f = [x_seg_f; x_f(i)];
        x_seg_b = [x_b(i); x_seg_b];
    end
    
    x_predict = a.*x_f + (1-a).*x_b;
    x_in(zero_location:zero_location+99) = x_predict';

end

figure; hold all;
plot(x_OG)
plot(x_in)