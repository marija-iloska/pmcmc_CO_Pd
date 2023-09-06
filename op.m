clear all
close all
clc

% Particle Metropolis within Gibbs
load area_ref490.mat
load expected_coverage.mat

% System specifications
tp_idx = 45;
cut_off = 0.30;
t_idx = 6;

% Data
time = time_mat_area{t_idx};
T = length(time);
y = area{t_idx};

% Some priors
eps_sat = mean(y(tp_idx - 25 : tp_idx))/cov_sat(t_idx);

% Number of particles
M = 50;


% Noise
var_A = (std(y(T-10:T)))^2;


% Bounds for state theta (coverage)
theta_max = 0.5;
theta_min = 0;


% System specifications
sys_specs = {var_A, eps_sat, cov_sat(t_idx)};
bounds = {tp_idx, cut_off, theta_max, theta_min};
P = 0.001;
dt = 0.067;

% METROPOLIS HASTINGS
alpha4 = 2;
beta4 = 2;

alpha3 = 3;
beta3 = 3;

alpha2 = 1000;
beta2 = 5;

alpha1 = 1000;
beta1 = 5;

% Propose k4
k4 = betarnd(alpha4, beta4);
k3 = betarnd(alpha3, beta3);


% Propose k1
x1 = gamrnd(alpha1, beta1);
k1 = k4*cut_off/(theta_max - cut_off)/P + x1;

% Propose k2
x2 = gamrnd(alpha1, beta1);
k2 = k3*cov_sat(t_idx)/(theta_max - cov_sat(t_idx))/P + x2;

% Prep parameters
a1 = k1*dt*P;
a4 = 1 - k4*dt;
a2 = k2*dt*P;
a3 = 1 - k3*dt;

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
J = 1000;
J0 = round(J/2);

tic
for j = 1:J


    % SAMPLE entire PF
    [theta_sample, epsilon_sample] = pf_chem(y, sys_specs, bounds, a, b, M, tp_AB, alpha);
    theta_chain(j,:) = theta_sample;
    epsilon_chain(j,:) = epsilon_sample;

  
    if (j > 3)
        tp_AB = find(theta_sample > cut_off);
        tp_AB = [tp_AB(1), tp_AB(end)];
        regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx + 1 : tp_AB(2), tp_AB(2):T};
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
    beta_A = 1./var_A + 0.5*sum( (y' - theta_sample.*epsilon_sample).^2 );
    var_a(j) = 1./gamrnd(1, beta_A);
    sys_specs = {var_a(j), eps_sat, cov_sat(t_idx)};
     

end
toc

% J0 = 3000;

theta_est = mean(theta_chain(J0:J,:),1);
epsilon_est = mean(epsilon_chain(J0:J,:),1);


figure;
plot(time, epsilon_est.*theta_est)
hold on
plot(time, y)
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
plot(x23chain(1:J,1), 'linewidth', 1)
title('Regions 2 and 3', 'FontSize', 15)

figure;
plot(x23chain(1:J,2), 'linewidth', 1)
title('Regions 2 and 3', 'FontSize', 15)


figure;
plot(x14chain(1:J,1), 'k', 'linewidth', 1)
title('Regions 1 and 4', 'FontSize', 15)

figure;
plot(x14chain(1:J,2), 'k', 'linewidth', 1)
title('Regions 1 and 4', 'FontSize', 15)

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
 
%save('Data/470op.mat')

