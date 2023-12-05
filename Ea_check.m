close all
clear all
clc

% RUN Ea
load Data/temps_info.mat

% temps_strings(4)=[];
N = length(temps_strings);

for n = 1 : N
    str = join(['Results/', temps_strings{n}, 'K.mat']);
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

J0 = 2000;
% Which data point to exclude
idx = setdiff(1:N, [1]);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;
% [Ea4, A4_des, Ea4_SE, A4_SE, ln_k4, Rsq_Ea] = get_Ea(k4_desorb(idx), T(idx), R);
% [Ea3, A3_des, Ea3_SE, A3_SE, ln_k3, Rsq_Ea] = get_Ea(k3_desorb(idx), T(idx), R);
% 
% [Ea2, A2_ads, Ea2_SE, A2_SE, ln_k2, Rsq_Ea] = get_Ea(k2_adsorb(idx), T(idx), R);
% [Ea1, A1_ads, Ea1_SE, A1_SE, ln_k1, Rsq_Ea] = get_Ea(k1_adsorb(idx), T(idx), R);

for j = 1:(J-J0-1)
    [Ea4(j), A4_des(j), Ea4_SE, A4_SE, ln_k4(:,j), Rsq_Ea] = get_Ea(kk4(idx,J0+j), T(idx), R);
    [Ea3(j), A3_des(j), Ea3_SE, A3_SE, ln_k3(:,j), Rsq_Ea] = get_Ea(kk3(idx,J0+j), T(idx), R);
    
    [Ea2(j), A2_ads(j), Ea2_SE, A2_SE, ln_k2(:,j), Rsq_Ea] = get_Ea(kk2(idx,J0+j), T(idx), R);
    [Ea1(j), A1_ads(j), Ea1_SE, A1_SE, ln_k1(:,j), Rsq_Ea] = get_Ea(kk1(idx,J0+j), T(idx), R);
end

Ea_mean = mean([Ea1; Ea2; Ea3; Ea4],2);
Ea = {Ea1, Ea2, Ea3, Ea4};
A = {A1_ads, A2_ads, A3_des, A4_des};

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
plot(1./T(idx), ln_kk1)

subplot(2,2,2)
scatter(1./T(idx), lnkk2(idx), 'filled')
hold on
plot(1./T(idx), ln_kk2)

subplot(2,2,3)
scatter(1./T(idx), lnkk3(idx), 'filled')
hold on
plot(1./T(idx), ln_kk3)

subplot(2,2,4)
scatter(1./T(idx), lnkk4(idx), 'filled')
hold on
plot(1./T(idx), ln_kk4)


% figure;
% subplot(2,2,1)
% scatter(1./T(idx), log(k1_adsorb(idx)), 'filled')
% hold on
% plot(1./T(idx), ln_k1)
% 
% subplot(2,2,2)
% scatter(1./T(idx), log(k2_adsorb(idx)), 'filled')
% hold on
% plot(1./T(idx), ln_k2)
% 
% subplot(2,2,3)
% scatter(1./T(idx), log(k3_desorb(idx)), 'filled')
% hold on
% plot(1./T(idx), ln_k3)
% 
% subplot(2,2,4)
% scatter(1./T(idx), log(k4_desorb(idx)), 'filled')
% hold on
% plot(1./T(idx), ln_k4)

% figure;
% for n = 1:N-1
%     plot(time_mat_area{n}, theta{n}, 'linewidth',1)
%     hold on
% end


save('J5000_Ea.mat', 'T', 'lnkk', 'ln_kk', "kk", 'R', 'Ea', 'A', 'Ea_mean' )
