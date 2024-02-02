clear all
close all
clc

% Particle Metropolis within Gibbs
load Data/area_ref490.mat
load Data/expected_coverage.mat
load Data/epsilons490.mat
load Data/temps_info.mat
load Data\colors.mat


% System specifications
tp_idx = 45;
cut_off = 0.33;
t_idx = 6;

% Data
time = time_area{t_idx};
y = area{t_idx}';
T = length(y);
str = temps_strings{t_idx};

% % Some priors
eps_sat = mean(y(tp_idx - 22 : tp_idx+2))/cov_sat(t_idx);
eps_exp = epsilon_exp(t_idx);


% Number of particles
M = 60;


% Noise
var_A = (std(y(T-5:T)))^2;

mu_E14 = eps_exp;
mu_E23 = eps_sat;
var_E14 = 0.5e-5;
var_E23 = 0.5e-6;

e14 = eps_exp;
e23 = eps_sat;



% Bounds for state theta (coverage)
theta_max = 0.5;
theta_min = 0;


% System specifications
sys_specs = {var_A, eps_sat, cov_sat(t_idx), eps_exp};
bounds = {tp_idx, cut_off, theta_max, theta_min};
P = 0.001;
dt = 0.067;

% METROPOLIS HASTINGS
alpha4 = 3;
beta4 = 3;

alpha3 = 3;
beta3 = 3;

alpha2 = 1000;
beta2 = 2;

alpha1 = 1000;
beta1 = 2;

% Propose k4
k4 = betarnd(alpha4, beta4);
k3 = betarnd(alpha3, beta3)/dt;
a4 = 1 - k4*dt;
a3 = 1 - k3*dt;

% Propose k1
x1 = gamrnd(alpha1, beta1);
k1 = k4*cut_off/(theta_max - cut_off)/P + x1;
k1_lim = (M - a4*cut_off)/(dt*P*(M - cut_off));

k1 = min([k1_lim, k1]);

% Propose k2
x2 = gamrnd(alpha1, beta1);
k2 = k3*cov_sat(t_idx)/(theta_max - cov_sat(t_idx))/P + x2;
k2_lim = (M - a3*cov_sat)/(dt*P*(M - cov_sat));

k2 = min([k2_lim, k2]);


% Prep parameters
a1 = k1*dt*P;
a2 = k2*dt*P;


% MH inputs
x14 = [k1, k4, x1];
x23 = [k2, k3, x2];


% PF Input
a = [a1, a2, 0, 0];
b = [a4, a3, a3, a4];

tp_AB = [5, 60];
regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx + 1 : tp_AB(2), tp_AB(2):T};

alpha = 5;
var = 0.01;

% Run GIBBS
J = 5000;
J0 = round(J/2);

tic
for j = 1:J


    % SAMPLE entire PF
    [theta_sample] = pf_chem_E(y, sys_specs, bounds, a, b, M, tp_AB, alpha);
    theta_chain(j,:) = theta_sample;

  
    if (j > 3)
        tp_AB = find(theta_sample > cut_off);
        %tp_AB = [5, 60];
        tp_AB = [tp_AB(1), tp_AB(end)];
        regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx : tp_AB(2), tp_AB(2):T};
    end

    % Sample Region 1 and 4 
    x = MH14op(theta_sample, regions, x14, bounds, P, dt, var, alpha4, alpha, alpha1, beta1);
    k1 = x(1);
    k4 = x(2);
    x1 = x(3);
    x14 = [k1, k4, x1];
    x14chain(j,:) = x14;


    % Sample Region 2 and 3
    x = MH23op(theta_sample, regions, x23,  bounds, P, dt, var, alpha3, alpha, alpha2, beta2, cov_sat(t_idx));
    k2 = x(1);
    k3 = x(2);
    x2 = x(3);
    x23 = [k2, k3, x2];
    x23chain(j,:) = x23;


    
    % Prep parameters
    a1 = k1*dt*P;
    a4 = 1 - k4*dt;
    a2 = k2*dt*P;
    a3 = 1 - k3*dt;

    % Concatenate for PF Input
    a = [a1, a2, 0, 0];
    b = [a4, a3, a3, a4];


    % Posterior for measurement variance
    e1_sample = e14*ones(1,length(regions{1}));
    e4_sample = e14*ones(1,length(regions{4}));
    e2_sample = e23*ones(1,length(regions{2}));
    e3_sample = e23*ones(1,length(regions{3}));

    epsilon_sample = [e1_sample, e2_sample, e3_sample, e4_sample(1:end-1)];
    epsilon_sample =epsilon_sample(1:length(theta_sample));
    beta_A = 1./var_A + 0.5*sum( (y' - theta_sample.*epsilon_sample').^2 );
    var_a(j) = 1./gamrnd(1, beta_A);

    % Sample epsilons
    e14 = sample_epsilon(y([regions{1}, regions{4}]), theta_sample([regions{1}, regions{4}])', var_a(j), var_E14, mu_E14);
    e23 = sample_epsilon(y([regions{2}, regions{3}]), theta_sample([regions{2}, regions{3}])', var_a(j), var_E23, mu_E23);
    
    eps_sat = e23;
    eps_exp = e14;

    epsilon_chain(j,:) = epsilon_sample;

    eps12(j,:) = [e14, e23];

    sys_specs = {var_a(j), eps_sat, cov_sat(t_idx), eps_exp};
     

end
toc


theta_est = mean(theta_chain(J0:J,:),1);
epsilon_est = mean(epsilon_chain(J0:J,:),1);


figure;
plot(time(1:T), epsilon_est.*theta_est)
hold on
plot(time(1:T), y)
title('Epsilon', 'FontSize', 15)

figure;
plot(theta_chain(1,:))
hold on
plot(theta_chain(J0,:))
hold on
plot(theta_chain(J,:))
hold on
plot(theta_est, 'k', 'linewidth',2)


figure;
plot(eps12(:,1), 'r', 'LineWidth',1)
hold on
plot(eps12(:,2), 'b', 'LineWidth',1)

figure;
subplot(2,2,1)
plot(x14chain(1:J,1), 'Color',col{1}, 'linewidth', 1)
title('Adsorption k_1', 'FontSize', 15)

subplot(2,2,2)
plot(x14chain(1:J,2), 'Color',col{1},'linewidth', 1)
title('Desorption k_4', 'FontSize', 15)

subplot(2,2,3)
plot(x23chain(1:J,1), 'k', 'linewidth', 1)
title('Adsorption k_2 ', 'FontSize', 15)

subplot(2,2,4)
plot(x23chain(1:J,2), 'k', 'linewidth', 1)
title('Desorption k_3 ', 'FontSize', 15)


k4_est = mean(x14chain(J0:J, 2),1);
k1_est = mean(x14chain(J0:J, 1),1);

k3_est = mean(x23chain(J0:J, 2),1);
k2_est = mean(x23chain(J0:J, 1),1);


figure;
subplot(2,2,1)
hist(x23chain(J0:J,1))
title('R2 Ads', 'FontSize', 15)


subplot(2,2,2)
hist(x23chain(J0:J,2))
title('R3 Des', 'FontSize', 15)


subplot(2,2,3)
hist(x14chain(J0:J,1))
title('R1 Ads', 'FontSize', 15)


subplot(2,2,4)
hist(x14chain(J0:J,2))
title('R4 Des', 'FontSize', 15)

sgtitle(str, 'FontSize', 15)

filename = join(['Results/noise_', str,'K_J5000.mat']);
save(filename)

