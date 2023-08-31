function [theta_sample, epsilon_sample] = pf_chem(y, sys_specs, bounds, a, b, M, tp_AB, alpha)


% Variances
[var_A, eps_sat, cov_sat] = sys_specs{:};
[tp_idx, cut_off, theta_max, theta_min] = bounds{:};

r = 1;

% Length of data
T = length(y);

% Initialize particles
theta_particles = beta_random(alpha, 2*0.1*ones(1,M))/2;
%theta_particles = pertrnd(theta_min, 0.1*ones(1,M), theta_max);
epsilon_particles = exprnd(0.5, 1, M);

% Get first estimates
theta_est(1) = mean(theta_particles);
epsilon_est(1) = mean(epsilon_particles);



for t = 2:T


%     % Which region are we in
     mean_eps = {0.5, eps_sat, eps_sat, 0.5};
     theta_mean = {a(r)*(0.5 - theta_particles) + b(r)*theta_particles, cov_sat*ones(1,M), a(r)*(0.5 - theta_particles) + b(r)*theta_particles, a(r)*(0.5 - theta_particles) + b(r)*theta_particles};

    % Propose particles
    theta_particles = beta_random(alpha, 2*theta_mean{r})/2;
    %theta_particles = pertrnd(theta_min-eps, theta_mean{r}, theta_max+eps);
    epsilon_particles = exprnd(mean_eps{r}, 1,M);

    if (isnan(theta_particles))
        disp('stop')
    end



    % Compute epsilon weights
    [w_cov, w_eps, theta_est(t), epsilon_particles, theta_particles] = compute_weights(y(t), epsilon_particles, theta_particles, var_A, M);

    % Store samples
    theta_store(t,:) = theta_particles;
    epsilon_store(t,:) = epsilon_particles;  

    % Estimate
    %theta_est(t) = mean(theta_particles);
    epsilon_est(t) = mean(epsilon_particles.*w_eps);
    idx_eps = datasample(1:M, M, 'Weights', w_eps);
    epsilon_particles = epsilon_particles(idx_eps);

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
%     if (t > tp_AB(1))
%        r = 2;
%        if (t > tp_idx)
%            r = 3;
%        elseif (t > tp_AB(2))
%            r = 4;
%        end
% 
%     end

    
end

% Take one sample of entire Time horizon
idx = datasample(1:M, 1, 'Weights', w_cov);
theta_sample = theta_store(:, idx);

idx_eps = datasample(1:M, 1, 'Weights', w_eps);
epsilon_sample = epsilon_store(:, idx_eps);

