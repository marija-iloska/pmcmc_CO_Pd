clear all
close all
clc


% DATA and SYSTEM ======================================================

% Load: Area time series, Saturation coverage, Prior estimated epsilons
% ... Labels for temperatures, Plotting colors
load Data/my_areas.mat
load Data/expected_coverage.mat
load Data/my_epsilons.mat
load Data/temps_info.mat
load Data/colors.mat

% Temperature indices
% T_idx = {1: 450K, 2: 460K, 3: 470K, 4: 475K, 5: 480K, 6: 490K}

% System specifications
tp_idx = 45;        % Time index when P is set to 0
cut_off = 0.27;     % Coverage around which we expect phase transition
t_idx = 6;          % Temperature index to choose


% Data Notation
time = time_area{t_idx};
y = area{t_idx}';


% Length of data and temperature label
T = length(y);
str = temps_strings{t_idx};


% Prior estimations: Epsilon at saturation and lower regions
eps_sat = mean(y(tp_idx - 10 : tp_idx+1))/cov_sat(t_idx);
eps_exp = epsilon_exp(t_idx);

% Obtain epsilons between low and high regions for all points
w = (y./ max(y));
epsilon = w.*eps_sat + (1-w).*eps_exp;


% Physical Constraints: Bounds for state theta (coverage)
theta_max = 0.5;
theta_min = 0;

% Pressure in Pascals and time step length
P = 0.001;
dt = 0.067;


% Prior Noise: Sample deviation of final points
sigma_A = (std(y(T-5:T)))/5;

% System specifications
sys_specs = {sigma_A, epsilon, cov_sat(t_idx)};
bounds = {tp_idx, cut_off, cov_sat(t_idx), theta_max, theta_min};


% METROPOLIS HASTINGS (MH) ==================================================

% Gamma Priors

% Region 1 and 4
alpha1 = 500;
beta1 = 20;

alpha4 = 2;
beta4 = 2;

% Regions 2 and 3
alpha2 = 500;
beta2 = 20;

alpha3 = 2;
beta3 = 2;


% Initialize samples of rate parameters k's

% Propose k4 and k3
k4 = betarnd(alpha4, beta4)/dt;
k3 = betarnd(alpha3, beta3)/dt;

% Prep model params
a4 = 1 - k4*dt;
a3 = 1 - k3*dt;

% Propose k1
x1 = gamrnd(alpha1, beta1);
k1 = k4*cut_off/(theta_max - cut_off)/P + x1;
k1_lim = (theta_max - a4*cut_off)/(dt*P*(theta_max - cut_off));
k1 = min([k1_lim, k1]);

% Propose k2
x2 = gamrnd(alpha1, beta1);
k2 = k3*cov_sat(t_idx)/(theta_max - cov_sat(t_idx))/P + x2;
k2_lim = (theta_max - a3*cov_sat)/(dt*P*(theta_max - cov_sat));
k2 = min([k2_lim, k2]);


% Prep model parameters
a1 = k1*dt*P;
a2 = k2*dt*P;


% Format as input in MH function
x14 = [k1, k4, x1];
x23 = [k2, k3, x2];


% Format as input in particle filter (PF)
a = [a1, a2, 0, 0];
b = [a4, a3, a3, a4];

% Initialize region index boundaries (they will vary)
tp_AB = [5, 60];
regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx + 1 : tp_AB(2), tp_AB(2):T};

% Prior param for sampling theta
alpha_theta = 5;

% SAMPLER SETTINGS

% Number of Gibbs iterations and burn-in
I = 200;
I0 = round(I/2);

% Number of particles
M = 40;

% Initialize arrays
theta_chain = zeros(I,T);
x14chain = zeros(I, 3);
x23chain = x14chain;
sigma_Achain = zeros(1,I);


%% PROPOSED PMCMC SAMPLER 

tic
for i = 1:I


    % SAMPLE entire PF --------------------------------------------------------
    [theta_sample] = pf_chem(y, sys_specs, bounds, a, b, M, alpha_theta);
    theta_chain(i,:) = theta_sample;

    % Update region indices
    tp_AB = find(theta_sample > cut_off);
    tp_AB = [tp_AB(1), tp_AB(end)];
    regions = {1 : tp_AB(1), tp_AB(1)+1 : tp_idx, tp_idx : tp_AB(2), tp_AB(2):T};

    % Sample Region 1 and 4 ------------------------------------------------------------------
    x = MH14(theta_sample, regions, x14, bounds, P, dt, alpha4, alpha_theta, alpha1, beta1);
    k1 = x(1);
    k4 = x(2);
    x1 = x(3);
    x14 = [k1, k4, x1];
    x14chain(i,:) = x14;


    % Sample Region 2 and 3 -----------------------------------------------------------------------------------
    x = MH23(theta_sample, regions, x23,  bounds, P, dt, alpha3, alpha_theta, alpha2, beta2);
    k2 = x(1);
    k3 = x(2);
    x2 = x(3);
    x23 = [k2, k3, x2];
    x23chain(i,:) = x23;

    
    % Prep parameters
    a1 = k1*dt*P;
    a4 = 1 - k4*dt;
    a2 = k2*dt*P;
    a3 = 1 - k3*dt;

    % Concatenate for PF Input
    a = [a1, a2, 0, 0];
    b = [a4, a3, a3, a4];


    % Sample noise sigma -------------------------------------------------
    beta_A = 1./sigma_A + 0.5*sum( (y - theta_sample.*epsilon).^2 );
    var_a = 1./gamrnd(1, beta_A);
    sigma_Achain(i) = var_a;

    % Update input parameters
    sys_specs{1} = var_a;
     

end
toc

% Get estimates
theta_est = mean(theta_chain(I0:I,:),1);

k4_est = mean(x14chain(I0:I, 2),1);
k1_est = mean(x14chain(I0:I, 1),1);

k3_est = mean(x23chain(I0:I, 2),1);
k2_est = mean(x23chain(I0:I, 1),1);


%% Visualize Data Fitting and Chain convergence

figure;
% Area vs Estimated Area
subplot(3,1,1)
plot(time(1:T), y, 'LineStyle','--', 'LineWidth',1)
hold on
plot(time(1:T), epsilon.*theta_est', 'LineStyle','-', 'LineWidth',1)
title('Epsilon * Theta', 'FontSize', 15)
legend('Data Area', 'Estimated Area', 'FontSize',12)


% Chains of the Coverage and Estimate
subplot(3,1,2)
plot(theta_chain(1,:))
hold on
plot(theta_chain(I0,:))
hold on
plot(theta_chain(I,:))
hold on
plot(theta_est, 'k', 'linewidth',2)
xlabel('Time', 'FontSize',12)
xlabel('Coverage', 'FontSize',12)
legend('1st Sample', 'I_0th Sample', 'Final Sample', 'Final Estimate', 'FontSize',12)


% Chain of the noise
subplot(3,1,3)
plot(sigma_Achain, 'k', 'linewidth',2)
xlabel('Iterations', 'FontSize',12)
ylabel('Noise', 'FontSize',12)

%% k PARAMETER CHAINS

figure;

% Region I and IV
subplot(2,2,1)
plot(x14chain(1:I,1), 'Color',col{1}, 'linewidth', 1)
title('Adsorption k_1', 'FontSize', 15)

subplot(2,2,2)
plot(x14chain(1:I,2), 'Color',col{1},'linewidth', 1)
title('Desorption k_4', 'FontSize', 15)
ylim([0,0.1])

% Region II and III
subplot(2,2,3)
plot(x23chain(1:I,1), 'k', 'linewidth', 1)
title('Adsorption k_2 ', 'FontSize', 15)

subplot(2,2,4)
plot(x23chain(1:I,2), 'k', 'linewidth', 1)
title('Desorption k_3 ', 'FontSize', 15)
ylim([0,1.5])


%% k PARAMETER HISTOGRAMS

figure;

% Region II and III
subplot(2,2,1)
hist(x23chain(I0:I,1))
title('R2 Ads', 'FontSize', 15)

subplot(2,2,2)
hist(x23chain(I0:I,2))
title('R3 Des', 'FontSize', 15)

% Region I and IV
subplot(2,2,3)
hist(x14chain(I0:I,1))
title('R1 Ads', 'FontSize', 15)

subplot(2,2,4)
hist(x14chain(I0:I,2))
title('R4 Des', 'FontSize', 15)

sgtitle(str, 'FontSize', 15)

%% SAVE RESULTS

% str_iter = num2str(I);
% filename = join(['RESULTS/pmcmc_', str,'K_I', str_iter, '.mat']);
% save(filename)

