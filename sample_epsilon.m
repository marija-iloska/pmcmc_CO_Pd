function [eps_sample] = sample_epsilon(A, theta, var_A, var_E, mu_E)

% Compute LIKELIHOOD
mu_L = sum(A.*theta)/sum(theta);
var_L = var_A/sum(theta.^2);

% Compute posterior mean
mu_post = (var_E * mu_L + var_L * mu_E)/(var_L + var_E);
var_post = var_L*var_E/(var_L + var_E);

% Sample epsilon
eps_sample = normrnd(mu_post, var_post);


end