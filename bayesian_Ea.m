close all
clear all
clc

% RUN Ea
load Data/temps_info.mat
load Data/colors.mat

% Number of Temperature datasets
N = length(temps_strings);

% Extract sampled parameters for each temperature
for n = 1 : N

    % Create filename and load data
    str = join(['Results/my_noise', temps_strings{n}, 'K_J5000.mat']);
    load(str)

    % Store Estimates
    k1_adsorb(n) = k1_est;
    k2_adsorb(n) = k2_est;

    k3_desorb(n) = k3_est;
    k4_desorb(n) = k4_est;


end

% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;

%% REGION IV
clc
alpha = 3;
beta_a = 100;
x = -1./(R*T);

N_samples = 1000;


beta_post = ((1 - beta_a*sum(x.*log(-log(k4_desorb))))/beta_a); 

for i = 1:N_samples
    Ea4(i) = gamrnd(alpha, beta_post);
end

figure;
subplot(2,2,4)
hist(Ea4)
title('Region 4')
mean(Ea4)



% REGION III
beta_post = ((1 - beta_a*sum(x.*log(-log(k3_desorb))))/beta_a); 

for i = 1:N_samples
    Ea3(i) = gamrnd(alpha, beta_post);
end

subplot(2,2,3)
hist(Ea3)
title('Region 3')
mean(Ea3)


alpha = alpha*0.1;


% REGION II
beta_post = ((1 - beta_a*sum(x.*log(log(k2_adsorb))))/beta_a); 

for i = 1:N_samples
    Ea2(i) = gamrnd(alpha, beta_post);
end

subplot(2,2,2)
hist(Ea2)
title('Region 2')
mean(Ea2)



% REGION I
beta_post = ((1 - beta_a*sum(x.*log(log(k1_adsorb))))/beta_a); 
for i = 1:N_samples
    Ea1(i) = gamrnd(alpha, beta_post);
end

subplot(2,2,1)
hist(Ea1)
title('Region 1')
mean(Ea1)