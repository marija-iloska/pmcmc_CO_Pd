function [w_cov,w_eps, theta_est, epsilon_particles, theta_particles] = compute_weights(y, epsilon_particles, theta_particles, var_A, M)


% Compute epsilon weights
ln_w_cov = -0.5*log(2*pi*var_A) - 0.5*(  (y - mean(epsilon_particles)*theta_particles).^2 )/var_A;

% Scale and normalize
w_cov = exp(ln_w_cov - max(ln_w_cov));
w_cov(isnan(w_cov)) = 0;
w_cov = w_cov./sum(w_cov);

if (isnan(w_cov))
    disp('stop')
end

idx_cov = datasample(1:M, M, 'Weights', w_cov);
theta_particles = theta_particles(idx_cov);
theta_est = mean(theta_particles);


% EPS weights
ln_w_eps = -0.5*log(2*pi*var_A) - 0.5*( (y - epsilon_particles.*theta_est).^2)/var_A;


% Sample one epsilon particle
w_eps = exp(ln_w_eps - max(ln_w_eps));
w_eps = w_eps./sum(w_eps);

% idx_eps = datasample(1:M, M, 'Weights', w_eps);
% epsilon_particles = epsilon_particles(idx_eps);
% epsilon_est = sum(epsilon_particles);


end