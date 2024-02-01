close all
clear all
clc

% RUN Ea
load Data/temps_info.mat
load Data/colors.mat

% temps_strings(4)=[];
N = length(temps_strings);

for n = 2 : N
    str = join(['Results/new_', temps_strings{n}, 'K_J1000.mat']);
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

J0 = 150;
clear Ea1 Ea2 Ea3 Ea4

% Which data point to exclude
idx = setdiff(1:N, [1,6]);

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
hold on
xline(0.5, 'Color', 'r', 'linewidth',3)
hold on
xline(1.5, 'Color', 'r', 'linewidth',3)
hold on
scatter(Ea_mean(1), 0, 70, 'g', 'filled')
title('E_1', 'FontSize', 15)

subplot(2,2,2)
hist(Ea2)
hold on
scatter(Ea_mean(2), 0, 70, 'g', 'filled')
title('E_2', 'FontSize', 15)

subplot(2,2,3)
hist(Ea3)
hold on
scatter(Ea_mean(3), 0, 70, 'g', 'filled')
title('E_3', 'FontSize', 15)

subplot(2,2,4)
hist(Ea4)
hold on
xline(24, 'Color', 'r', 'linewidth',3)
hold on
xline(36, 'Color', 'r', 'linewidth',3)
scatter(Ea_mean(4), 0, 70, 'g', 'filled')
title('E_4', 'FontSize', 15)




% A   PREEXP FACTOR
figure;
subplot(2,2,1)
hist(ln_A1)
hold on
scatter(lnA_mean(1), 0, 70, 'g', 'filled')
title('ln(A_1)', 'FontSize', 15)

subplot(2,2,2)
hist(ln_A2)
hold on
scatter(lnA_mean(2), 0, 70, 'g', 'filled')
title('ln(A_2)', 'FontSize', 15)

subplot(2,2,3)
hist(ln_A3)
hold on
scatter(lnA_mean(3), 0, 70, 'g', 'filled')
title('ln(A_3)', 'FontSize', 15)

subplot(2,2,4)
hist(ln_A4)
hold on
xline(log(10^13.5), 'Color', 'r', 'linewidth',3)
scatter(lnA_mean(4), 0, 70, 'g', 'filled')
title('ln(A_4)', 'FontSize', 15)






figure;
subplot(2,2,1)
scatter(1./T(idx), lnkk1(idx), 'k','filled')
hold on 
plot(1./T(idx), ln_kk1(idx),'Color',  col{4}, 'linewidth',2)
title('Adsorption k_1', 'FontSize', 15)
ylim([7,9])
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

subplot(2,2,2)
scatter(1./T(idx), lnkk2(idx), 'k', 'filled')
hold on
plot(1./T(idx), ln_kk2(idx), 'Color', col{4}, 'LineWidth',2)
title('Adsorption k_2', 'FontSize', 15)
ylim([7,9])
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

subplot(2,2,3)
scatter(1./T(idx), lnkk3(idx), 'k', 'filled')
hold on
plot(1./T(idx), ln_kk3(idx), 'Color', col{1}, 'LineWidth', 2)
title('Desorption k_3', 'FontSize', 15)
ylim([-8,2])
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

subplot(2,2,4)
scatter(1./T(idx), lnkk4(idx),  'k','filled')
hold on
plot(1./T(idx), ln_kk4(idx), 'Color', col{1}, 'LineWidth', 2)
title('Desorption k_4', 'FontSize', 15)
ylim([-8,2])
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on






figure;
for n = 2:N
    plot(time_mat_area{n}(1:length(theta{n})), theta{n}, 'linewidth',1.5, 'Color', col{n-1})
    hold on
end
xlabel('Time [s]','FontSize', 15)
ylabel('Covereage [ML]', 'FontSize', 15)
legend(temps_strings{2:end}, 'FontSize', 15)
title('Inferred Latent States', 'FontSize', 15)
grid on


%save('eps_J10000_Ea.mat', 'T', 'lnkk', 'ln_kk', "kk", 'R', 'Ea', 'ln_A', 'Ea_mean', "lnA_mean" )
