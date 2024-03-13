function [theta_sample] = pf_chem(y, sys_specs, bounds, a, b, M, alpha)


% Extract settings and constraints
[var_A, epsilon, cov_sat] = sys_specs{:};
[tp_idx, cut_off, ~, ~, ~] = bounds{:};

% Start with Region I
r = 1;

% Length of data
T = length(y);

% Initialize particles
theta_particles = beta_random(alpha, 2*0.1*ones(1,M))/2;
theta_store = zeros(T,M);

% Get first estimates
theta_est = zeros(T,1);
theta_est(1) = mean(theta_particles);


% Start Filter
for t = 2:T

    % Which region are we in
    temp_mean = a(r)*(0.5 - theta_particles) + b(r)*theta_particles;

    % Constrain from overflow
    mean_min = min(temp_mean, cov_sat*ones(1,M));
    theta_mean = {temp_mean, mean_min, mean_min, temp_mean};

    % Propose particles
    theta_particles = beta_random(alpha, 2*theta_mean{r})/2;

    if (isnan(theta_particles))
        disp('Theta was sampled outside of bounds. Try again.')
    end

    % Compute epsilon weights
    [w_cov, theta_est(t), theta_particles] = compute_weights(y(t), epsilon(t), theta_particles, var_A, M);

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

