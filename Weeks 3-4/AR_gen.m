clear all; clc; close all;

% Create AR(4) Model
N=1000;
P=4;

% Define pole method for AR
pole(1)=0.99*exp(j*0.1*pi);
pole(3)=conj(pole(1));
pole(2)=1.005*exp(j*0.4*pi);
pole(4)=conj(pole(2));
a = poly(pole); %Makes polynomial with those poles to produce the AR model

% Put actual AR PArameters into vector
theta = -a(2:end);

% Put actual AR PArameters into vector
theta = -a(2:end);

% Error
m_theta = 0; var_e=1;
e = m_theta + randn(N,1)*var_e;

% Past data points assumed to be zero:
x=filter(1,a,e);

% Plot data
figure
subplot(211), plot(e)
title('e_n')
subplot(212), plot(x)
title('x_n')

% Make a Matlab system model:
figure;
subplot(311), plot(x)
title('x_n')
sys=tf(1,a,1);
subplot(313), pzplot(sys)
title('Poles of AR(4) model')
[H,w]=freqz(1,a);
subplot(312),
semilogy(w,abs(H).^2*var_e^2)
title('Power spectrum of AR(4) model')
xlabel('Normalised Frequency (rads^{-1})')


% Form linear model:
g1 = [0 x(1:N-1)']';
g2 = [0 0 x(1:N-2)']';
g3 = [0 0 0 x(1:N-3)']';
g4 = [0 0 0 0 x(1:N-4)']';
G  = [g1 g2 g3 g4];

% Compute ML estimate
theta_ML = inv(G'*G)*G'*x;

% Compute MAP estimate
m_theta = [0; 0; 0; 0]; cov_theta = eye(4);
theta_MAP = inv(G'*G + var_e*inv(cov_theta))*(G'*x + var_e*inv(cov_theta)*m_theta);

% Compute MSE
x_ML  = filter(1,[1 -theta_ML'], e);
x_MAP = filter(1,[1 -theta_MAP'],e);

e_ML = sum((x_ML-x).^2)/N;
e_MAP = sum((x_MAP-x).^2)/N;
