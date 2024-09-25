function [theta_est, theta_particles] = compute_weights(y, epsilon, theta_particles, var_A, M)


% Compute epsilon weights
ln_w_cov = -0.5*log(2*pi*var_A) - 0.5*(  (y - epsilon*theta_particles).^2 )/var_A;

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


end