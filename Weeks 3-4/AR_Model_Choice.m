clear all; clc; close all;

% Create AR(4) Model
N=1000;
P=4;

% Define pole method for AR
pole(1)=0.99*exp(j*0.5*pi);
pole(3)=conj(pole(1));
pole(2)=0.97*exp(j*0.4*pi);
pole(4)=conj(pole(2));
a = poly(pole); %Makes polynomial with those poles to produce the AR model

% Put actual AR PArameters into vector
theta = -a(2:end);

% Put actual AR PArameters into vector
theta = -a(2:end);

% Error
m_e = 0; var_e=1;
e = m_e + randn(N,1)*var_e;

% Past data points assumed to be zero:
x=filter(1,a,e);
y = x;

% Examine a choice of 4 possible models (4 orders)
% AR(1)
m_1 = 0; cov_1 = 1;
G1 = [0 x(1:N-1)']';

% AR(2)
m_2 = 0; cov_2 = eye(2);
g1 = [0 x(1:N-1)']'; g2 = [0 0 x(1:N-2)']';
G2  = [g1 g2];

% AR(3)
m_3 = 0; cov_3 = eye(3);
g1 = [0 x(1:N-1)']'; g2 = [0 0 x(1:N-2)']'; g3 = [0 0 0 x(1:N-3)']';
G3  = [g1 g2 g3];

% AR(4)
m_4 = 0; cov_4 = eye(4);
g1 = [0 x(1:N-1)']'; g2 = [0 0 x(1:N-2)']'; g3 = [0 0 0 x(1:N-3)']'; g4 = [0 0 0 0 x(1:N-4)']';
G4  = [g1 g2 g3 g4];

% Generate MAP matrices
phi1 = G1'*G1 + var_e*inv(cov_1);
The1 = G1'*y  + var_e*inv(cov_1)*m_1;
phi2 = G2'*G2 + var_e*inv(cov_2);
The2 = G2'*y  + var_e*inv(cov_2)*m_2;
phi3 = G3'*G3 + var_e*inv(cov_3);
The3 = G3'*y  + var_e*inv(cov_3)*m_3;
phi4 = G4'*G4 + var_e*inv(cov_4);
The4 = G4'*y  + var_e*inv(cov_4)*m_3;

% MAP estimates of parameters
theta1_MAP = inv(phi1)*The1;
theta2_MAP = inv(phi2)*The2;
theta3_MAP = inv(phi3)*The3;
theta4_MAP = inv(phi4)*The4;

% Calulate each of the terms for the ML
P1 = 1; P2 = 2; P3 = 3; P4 = 4;

front1 = 1/((2*pi)^(P1/2) * sqrt(det(cov_1)) * sqrt(norm(phi1)) * (2*pi*var_e)^((N-P1)/2));
front2 = 1/((2*pi)^(P2/2) * sqrt(det(cov_2)) * sqrt(norm(phi2)) * (2*pi*var_e)^((N-P2)/2));
front3 = 1/((2*pi)^(P3/2) * sqrt(det(cov_3)) * sqrt(norm(phi3)) * (2*pi*var_e)^((N-P3)/2));
front4 = 1/((2*pi)^(P4/2) * sqrt(det(cov_4)) * sqrt(norm(phi4)) * (2*pi*var_e)^((N-P4)/2));

exp1 = -(1/(2*var_e)) * (y'*y + var_e * m_1' * inv(cov_1) * m_1 - The1'*theta1_MAP);
exp2 = -(1/(2*var_e)) * (y'*y + var_e * m_2' * inv(cov_2) * m_2 - The2'*theta2_MAP);
exp3 = -(1/(2*var_e)) * (y'*y + var_e * m_3' * inv(cov_3) * m_3 - The3'*theta3_MAP);
exp4 = -(1/(2*var_e)) * (y'*y + var_e * m_3' * inv(cov_4) * m_3 - The4'*theta4_MAP);

pM1 = front1 * exp(exp1); pM1 = pM1(1,1);
pM2 = front2 * exp(exp2); pM2 = pM2(1,1);
pM3 = front3 * exp(exp3); pM3 = pM3(1,1);
pM4 = front4 * exp(exp4); pM4 = pM4(1,1);

p_norm = pM1 + pM2 + pM3 + pM4;
pM1 = pM1/p_norm
pM2 = pM2/p_norm
pM3 = pM3/p_norm
pM4 = pM4/p_norm