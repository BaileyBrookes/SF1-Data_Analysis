clear all; clc; close all;

% Create AR(4) Model
N=100;
P_order=4;

% Define pole method for AR
pole(1)=0.99*exp(j*0.5*pi);
pole(3)=conj(pole(1));
pole(2)=0.97*exp(j*0.1*pi);
pole(4)=conj(pole(2));
a = poly(pole); %Makes polynomial with those poles to produce the AR model

% Put actual AR Parameters into vector
theta = -a(2:end);


% [x_in, Fs] = audioread('Sound Files/armst_37_orig.wav.'); 
% x = x_in;
% N = length(x);
% % Error
m_e = 0; var_e=5;
e = m_e + randn(N,1)*var_e;
% x = x ;

% Past data points assumed to be zero:
x=filter(1,a,e);

figure; hold on
for P = 1:8
    
    for i = 1:P
        G(:,i) = [zeros(i,1); x(1:N-i)];
    end
    
    m_theta = zeros(P,1); cov_theta = 5*eye(P); 
    phi = G'*G + var_e*inv(cov_theta);
    The = G'*x  + var_e*inv(cov_theta)*m_theta;
    theta_MAP{P} = inv(phi)*The;
    ML{P} = inv(G'*G)*G'*x;
    
    front = 1/((2*pi)^(P/2) * sqrt(det(cov_theta)) * sqrt(det(phi)) * (2*pi*var_e)^((N-P)/2));
    exp_term = -(1/(2*var_e)) * (x'*x + var_e * m_theta' * inv(cov_theta) * m_theta - The'*theta_MAP{P});
    pM = front * exp(exp_term); pM(P) = pM(1,1);
    
    % Compute MSE
    x_MAP{P} = G * theta_MAP{P};
    e_MAP(P) = sum((x_MAP{P}-x).^2)/N;
    
    a = [1 -theta_MAP{P}'];
    [H,w]=freqz(1,a);
    semilogy(w,abs(H).^2)
    title(['Power spectrum of AR(', num2str(P), ') model'])
end
figure; plot(e_MAP);
xlabel('AR Model Order P')
ylabel('Mean Squared Error')