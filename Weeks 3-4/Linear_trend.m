% Code for week 3 and 4 question 2

clear all; clc; close all;

% Generator Matrices
N = 100;
n = [1:N]';
G = [ones(N,1), n];

% Error
m_e =0; var_e = 1;
error = m_e + var_e*randn(1,N);

% Theta
m_theta = [5;5];
cov_theta= [100 0; 0 100];
theta = [1;2];

% Generate data
y = G*theta + error';

% Compute ML estimate
theta_ML = inv(G'*G)*G'*y;

% Comoute MAP esimate
theta_MAP = inv(G'*G + var_e*inv(cov_theta))*(G'*y + var_e*inv(cov_theta)*m_theta);