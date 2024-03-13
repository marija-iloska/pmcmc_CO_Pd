function [theta_sample] = pf_chem(y, sys_specs, bounds, a, b, M, tp_AB, alpha)


% Variances
%[var_A, eps_sat, cov_sat, eps_exp] = sys_specs{:};
[var_A, epsilon, cov_sat] = sys_specs{:};

[tp_idx, cut_off, ~, ~] = bounds{:};

r = 1;

% Length of data
T = length(y);

% Initialize particles
theta_particles = beta_random(alpha, 2*0.1*ones(1,M))/2;

% Get first estimates
theta_est(1) = mean(theta_particles);


for t = 2:T

     % Which region are we in
    mean_eps = epsilon(t);
    temp_mean = a(r)*(0.5 - theta_particles) + b(r)*theta_particles; 
    mean_min = min(temp_mean, cov_sat*ones(1,M));
    theta_mean = {temp_mean, mean_min, mean_min, temp_mean};

    % Propose particles
    theta_particles = beta_random(alpha, 2*theta_mean{r})/2;

    if (isnan(theta_particles))
        disp('stop')
    end


    % Compute epsilon weights
    [w_cov, theta_est(t), theta_particles] = compute_weights(y(t), mean_eps, theta_particles, var_A, M);

    % Store samples
    theta_store(t,:) = theta_particles;

    % Identify region
    if (theta_est(t) > cut_off)
        r = 2;
        if (t > tp_idx)
            r = 3;
        end
    else
        if (t > tp_idx)
            r = 4;
        end
    end

    
end

% Take one sample of entire Time horizon
idx = datasample(1:M, 1, 'Weights', w_cov);
theta_sample = theta_store(:, idx);

