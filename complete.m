clear all
close all
clc

% Particle Metropolis within Gibbs
load area_ref490.mat
load expected_coverage.mat

% System specifications
tp_idx = 45;
cut_off = 0.33;

t_idx = 1;

% Data
time = time_mat_area{t_idx};
T = length(time);
y = area{t_idx};

% Some priors
eps_sat = mean(y(tp_idx - 30 : tp_idx))/cov_sat(t_idx);

% Number of particles
M = 100;


% Noise
var_A = (std(y(T-10:T)))^2;


% Bounds for state theta (coverage)
theta_max = 0.5;
theta_min = 0;


% System specifications
sys_specs = {var_A, eps_sat, cov_sat(t_idx)};
bounds = {tp_idx, cut_off, theta_max, theta_min};

% METROPOLIS HASTINGS
alpha1 = 2;
beta1 = 2;

alpha2 = 2;
beta2 = 2;

alpha3 = 2;
beta3 = 2;

alpha4 = 2;
beta4 = 2;

% Propose sample for a
a4 = betarnd(alpha4, beta4);
y1 = betarnd(alpha1, beta1);

% Max data point in R1
y1max = 0.33;
y2max = cov_sat(t_idx);

% Upper bound of a1
a1_lim = (cut_off - a4*y1max)/(M - y1max);

% Sample a1
a1 = a1_lim*y1;

% Propose sample for a
a3 = betarnd(alpha3, beta3);
y2 = betarnd(alpha2, beta2);

% Upper bound of a1
a2_lim = (M - a4*y2max)/(M - y2max);

% Sample a1
a2 = a2_lim*y2;

% Concatenate
a14 = [a1, a4];
a23 = [a2, a3];

% Input
a = [a1, a2, 0, 0];
b = [a4, a3, a3, a4];

tp_AB = [5, 60];
regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx + 1 : tp_AB(2), tp_AB(2):T};

alpha = 5;

% Run GIBBS
J = 600;
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
    x = MH14(theta_sample, regions, a14, bounds, alpha1, alpha4, beta1, alpha);
    a14 = x;
    achain14(j,:) = a14;


    % Sample Region 2 and 3
    x = MH23(theta_sample, regions, a23,  bounds, alpha2, alpha3, beta2, alpha, cov_sat(t_idx));
    a23 = x;
    achain23(j,:) = a23;


    % Concatenate
    a = [a14(1), a23(1), 0, 0];
    b = [a14(2), a23(2), a23(2), a14(2)];


    % Posterior for measurement variance
    beta_A = 1./var_A + 0.5*sum( (y' - theta_sample.*epsilon_sample).^2 );
    var_a(j) = 1./gamrnd(1, beta_A);
    sys_specs = {var_a(j), eps_sat, cov_sat(t_idx)};
     

end
toc

theta_est = mean(theta_chain(J0:J,:),1);
epsilon_est = mean(epsilon_chain(J0:J,:),1);
% k_des = mean(kXO(J0:J,:),1);
% k_ads = mean(kOX(J0:J,:),1);

k_des = mean((1 - achain14(J0:J, 2))/0.067,1);

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
plot(achain23(:,1), 'linewidth', 1)

figure;
plot(achain23(:,2), 'k', 'linewidth', 1)
title('Regions 2 and 3', 'FontSize', 15)


figure;
plot(achain14(:,1), 'linewidth', 1)

figure;
plot(achain14(:,2), 'k', 'linewidth', 1)
title('Regions 1 and 4', 'FontSize', 15)

% figure;
% subplot(2,2,1)
% hist(kOX(:,1))
% title('R1 Ads', 'FontSize', 15)
% 
% 
% subplot(2,2,2)
% hist(kOX(:,2))
% title('R2 Ads', 'FontSize', 15)
% 
% 
% subplot(2,2,3)
% hist(kXO(:,1))
% title('R3 Des', 'FontSize', 15)
% 
% 
% subplot(2,2,4)
% hist(kXO(:,2))
% title('R4 Des', 'FontSize', 15)
% 
% save('450test.mat', 'k_des', 'k_ads')

