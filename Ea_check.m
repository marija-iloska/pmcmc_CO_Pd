close all
clear all
clc

% RUN Ea
load Data/temps_info.mat

% temps_strings(4)=[];
N = length(temps_strings);

for n = 1 : N
    str = join(['Results/eps_', temps_strings{n}, 'K_J10000.mat']);
    load(str)

    k1_adsorb(n) = k1_est;
    k2_adsorb(n) = k2_est;

    k3_desorb(n) = k3_est;
    k4_desorb(n) = k4_est;

    theta{n} = theta_est;

    kk1(n,:) = x14chain(:,1);
    kk2(n,:) = x23chain(:,1);
    kk3(n,:) = x23chain(:,2);
    kk4(n,:) = x14chain(:,2);

end




%% Plot

% Temperatures
T = [450, 460, 470, 475, 480, 490];

J0 = 7000;
clear Ea1 Ea2 Ea3 Ea4

% Which data point to exclude
idx = setdiff(1:N, [1]);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;


for j = 1:(J-J0-1)
    [Ea4(j), ln_A4(j), ln_k4(idx,j)] = get_Ea(kk4(idx,J0+j), T(idx), R);
    [Ea3(j), ln_A3(j),  ln_k3(idx,j)] = get_Ea(kk3(idx,J0+j), T(idx), R);
    
    [Ea2(j), ln_A2(j),  ln_k2(idx,j)] = get_Ea(kk2(idx,J0+j), T(idx), R);
    [Ea1(j), ln_A1(j),  ln_k1(idx,j)] = get_Ea(kk1(idx,J0+j), T(idx), R);
end

Ea_mean = mean([Ea1; Ea2; Ea3; Ea4],2)
lnA_mean = mean([ln_A1; ln_A2; ln_A3; ln_A4],2)
Ea = {Ea1, Ea2, Ea3, Ea4};
ln_A = {ln_A1, ln_A2, ln_A3, ln_A4};

ln_kk1 = mean(ln_k1, 2);
ln_kk2 = mean(ln_k2, 2);
ln_kk3 = mean(ln_k3, 2);
ln_kk4 = mean(ln_k4, 2);

lnkk1 = mean(log(kk1), 2);
lnkk2 = mean(log(kk2), 2);
lnkk3 = mean(log(kk3), 2);
lnkk4 = mean(log(kk4), 2);

kk = {kk1, kk2, kk3, kk4};
ln_kk = {ln_kk1, ln_kk2, ln_kk3, ln_kk4};
lnkk = {lnkk1, lnkk2, lnkk3, lnkk4};



figure;
subplot(2,2,1)
hist(Ea1)
title('Ea1')

subplot(2,2,2)
hist(Ea2)
title('Ea2')

subplot(2,2,3)
hist(Ea3)
title('Ea3')

subplot(2,2,4)
hist(Ea4)
title('Ea4')



figure;
subplot(2,2,1)
scatter(1./T(idx), lnkk1(idx), 'filled')
hold on
plot(1./T(idx), ln_kk1(idx))

subplot(2,2,2)
scatter(1./T(idx), lnkk2(idx), 'filled')
hold on
plot(1./T(idx), ln_kk2(idx))

subplot(2,2,3)
scatter(1./T(idx), lnkk3(idx), 'filled')
hold on
plot(1./T(idx), ln_kk3(idx))

subplot(2,2,4)
scatter(1./T(idx), lnkk4(idx), 'filled')
hold on
plot(1./T(idx), ln_kk4(idx))




figure;
for n = 2:N
    plot(time_mat_area{n}(1:length(theta{n})), theta{n}, 'linewidth',1)
    hold on
end


%save('eps_J10000_Ea.mat', 'T', 'lnkk', 'ln_kk', "kk", 'R', 'Ea', 'ln_A', 'Ea_mean', "lnA_mean" )
