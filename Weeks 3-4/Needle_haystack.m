% Bayesian Model Choice
clear all; close all; clc;

% Load data
load('hidden_data.mat')
N = length(y);

% Model 1: y = e
G1 = zeros(N,1);
theta1 = 0;
m_theta1 = 0; cov_theta1 = 0.01;
m_e1 = 0; var_e1 = 1; 

% Model 2: y = a + e
G2 = ones(N,1);
theta2 = 1;
m_theta2 = 1; cov_theta2 = 1;
m_e2 = 0; var_e2 = 1;

% Model 3: y a + bn + e
n = [1:N]';
G3 = [ones(N,1), n];
theta3 = [1;2];
m_theta3 = [1;2]; cov_theta3= [1 0; 0 1];
m_e3 =0; var_e3 = 1;

% Generate MAP matrices
phi1 = G1'*G1 + var_e1*inv(cov_theta1);
The1 = G1'*y  + var_e1*inv(cov_theta1)*m_theta1;
phi2 = G2'*G2 + var_e2*inv(cov_theta2);
The2 = G2'*y  + var_e2*inv(cov_theta2)*m_theta2;
phi3 = G3'*G3 + var_e3*inv(cov_theta3);
The3 = G3'*y  + var_e3*inv(cov_theta3)*m_theta3;

% MAP estimates of parameters
theta1_MAP = inv(phi1)*The1;
theta2_MAP = inv(phi2)*The2;
theta3_MAP = inv(phi3)*The3;