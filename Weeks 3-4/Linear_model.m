% Code for week 3 and 4 question 1

clear all; close all; clc;

% Generate data
N = 5;
theta = 4;
m_e = 0; var_e = 1;
error = m_e + var_e*randn(N,1);
G = ones(N,1);
y = G*theta + error;

% Compute ML estimate of theta
theta_ML = sum(y)/N;

% Compute the MAP estimate of theta
m_theta = 2; cov_theta = 0.5; % Prior on theta is a gaussian with these parameters
theta_MAP = inv(G'*G + var_e*inv(cov_theta))*(G'*y + var_e*inv(cov_theta)*m_theta);

% the maths way
var_frac = var_e/cov_theta;
theta_MAP2 = (sum(y)/N + m_theta*var_frac)/(1 + var_frac);

% Densities
prior     = makedist('Normal','mu',m_theta,    'sigma',cov_theta);
liklihood = makedist('Normal','mu',theta_ML,      'sigma',var_e);
posterior = makedist('Normal','mu',theta_MAP,  'sigma',var_e*inv((G'*G + var_e*inv(cov_theta))));

% Make pdfs
x = -5:0.01:15;
prior_pdf     = pdf(prior,    x);
liklihood_pdf = pdf(liklihood,x); 
posterior_pdf = pdf(posterior,x);

% Plot pdfs
figure; hold on;
plot(x,prior_pdf);
plot(x,liklihood_pdf);
plot(x,posterior_pdf);
xlabel('\theta')
ylabel('Probability Density')
legend('Prior: p(\theta)', 'Likelihood: p(y|\theta)', 'Posterior: p(\theta|y)','Location', 'NorthWest')