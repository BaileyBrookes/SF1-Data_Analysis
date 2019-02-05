clear all; close all; clc;

% Read data
[x_in, Fs] = audioread('Sound Files/grosse_original.wav'); 
x = x_in;
x(5780:5780+100) = 0;
N = length(x);

% Error
m_e = 0; var_e=0;

figure; hold on
for P = 1:100
    
    % Form generator matrix
    for i = 1:P
        G(:,i) = [zeros(i,1); x(1:N-i)];
    end
    
    % Form parpameter estimate fro MAP and ML
    m_theta = zeros(P,1); cov_theta = 5*eye(P);
    phi = G'*G + var_e*inv(cov_theta);
    The = G'*x  + var_e*inv(cov_theta)*m_theta;
    theta_MAP{P} = inv(phi)*The;
    theta_ML{P} = inv(G'*G)*G'*x;
    
    % Compute likihood of model
    front = 1/((2*pi)^(P/2) * sqrt(det(cov_theta)) * sqrt(det(phi)) * (2*pi*var_e)^((N-P)/2));
    exp_term = -(1/(2*var_e)) * (x'*x + var_e * m_theta' * inv(cov_theta) * m_theta - The'*theta_MAP{P});
    pM = front * exp(exp_term); pM(P) = pM(1,1);
    
    % Compute MSE
    x_MAP{P} = G * theta_MAP{P};
    e_MAP(P) = sum((x_MAP{P}-x).^2)/N;
    
    % Plot PSD for each model
    a = [1 -theta_MAP{P}'];
    [H,w]=freqz(1,a);
    semilogy(w,abs(H).^2)
    title(['Power spectrum of AR(', num2str(P), ') model'])
end
% Plot MSE vs Number of Parameters
figure; plot(e_MAP);

% Now choose a model
theta = theta_MAP{15}
x_MAP = x_MAP{15}
figure; hold on;
plot(x_MAP)
plot(x)

