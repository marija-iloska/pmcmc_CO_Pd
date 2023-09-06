close all
clear all
clc

% RUN Ea
load temps_info.mat

% temps_strings(4)=[];
N = length(temps_strings);

for n = 1 : N
    str = join(['Data/', temps_strings{n}, 'op.mat']);
    load(str)
    k1_adsorb(n) = k1_est;
    k2_adsorb(n) = k2_est;

    k3_desorb(n) = k3_est;
    k4_desorb(n) = k4_est;

    theta{n} = theta_est;

end


% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Which data point to exclude
idx = setdiff(1:N, [4]);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;


[Ea4_des, A4_des, Ea4_SE, A4_SE, ln_k, Rsq_Ea] = get_Ea(k4_desorb(idx), T(idx), R);
[Ea3_des, A3_des, Ea3_SE, A3_SE, ln_k, Rsq_Ea] = get_Ea(k3_desorb(idx), T(idx), R);

[Ea2_ads, A2_ads, Ea2_SE, A2_SE, ln_k, Rsq_Ea] = get_Ea(k2_adsorb(idx), T(idx), R);
[Ea1_ads, A1_ads, Ea1_SE, A1_SE, ln_k, Rsq_Ea] = get_Ea(k1_adsorb(idx), T(idx), R);

figure;
scatter(1./T(idx), log(k4_desorb(idx)), 'filled')

figure;
for n = 1:N-1
    plot(time_mat_area{n}, theta{n}, 'linewidth',1)
    hold on
end

