function [Y] = beta_random(alpha, mu)

% Compute beta parameter based on define alpha and expected mean
beta  = (1-mu).*alpha./mu;

% Sample random variable
Y = betarnd(alpha, beta);


end