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
    %str = join(['RESULTS/pmcmc_', temps_strings{n},'K_I3000', '.mat']);
    load(str)

    % Store Estimates
    k1_adsorb(n) = k1_est;
    k2_adsorb(n) = k2_est;

    k3_desorb(n) = k3_est;
    k4_desorb(n) = k4_est;

    theta{n} = theta_est;

    % Store chains
    kk1(n,:) = x14chain(:,1);
    kk2(n,:) = x23chain(:,1);
    kk3(n,:) = x23chain(:,2);
    kk4(n,:) = x14chain(:,2);

    % Only For plotting
    if n==2
        k2chain = x23chain(:,1);
        k4chain = x14chain(:,2);  
    end

    sig_store{n} = var_a;

end

I = J;


%% PLOTS

% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Burn-in
I0 = 2000;

% In case we want to exclude a tempereature
idx = setdiff(1:N, []);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;

% Compute Ea samples for each region using k samples (after burn-in)
for j = 1:(I-I0-1)
    [Ea4(j), ln_A4(j), ln_k4(idx,j)] = compute_Ea(kk4(idx,I0+j), T(idx), R);
    [Ea3(j), ln_A3(j),  ln_k3(idx,j)] = compute_Ea(kk3(idx,I0+j), T(idx), R);
    
    [Ea2(j), ln_A2(j),  ln_k2(idx,j)] = compute_Ea(kk2(idx,I0+j), T(idx), R);
    [Ea1(j), ln_A1(j),  ln_k1(idx,j)] = compute_Ea(kk1(idx,I0+j), T(idx), R);
end

% Compute mean estimates
Ea_mean = mean([Ea1; Ea2; Ea3; Ea4],2)
lnA_mean = mean([ln_A1; ln_A2; ln_A3; ln_A4],2)

Ea_freq = {Ea1, Ea2, Ea3, Ea4};
lnA_freq = {ln_A1, ln_A2, ln_A3, ln_A4};

%save('RESULTS/freq.mat', 'Ea_freq', "lnA_freq");


% Log fitting and estimates for plotting
ln_kk1 = mean(ln_k1, 2);
ln_kk2 = mean(ln_k2, 2);
ln_kk3 = mean(ln_k3, 2);
ln_kk4 = mean(ln_k4, 2);

lnkk1 = mean(log(kk1), 2);
lnkk2 = mean(log(kk2), 2);
lnkk3 = mean(log(kk3), 2);
lnkk4 = mean(log(kk4), 2);


% Store chains
% Ea = {Ea1, Ea2, Ea3, Ea4};
% ln_A = {ln_A1, ln_A2, ln_A3, ln_A4};
% ln_kk = {ln_kk1, ln_kk2, ln_kk3, ln_kk4};


% PLOTS ==================================================================

%% ACTIVATION ENERGY HISTOGRAMS 
% figure;
% for n = 1:4
%     p = subplot(2,2,n)
%     h = histogram(Ea_store{n});
%     h.FaceColor = [0, 0.35, 0.65];
%     h.EdgeColor = 'k'; % [0.8, 0.8, 0.8];
%     h.LineWidth = 0.01;
%     h.FaceAlpha = 0.95;
%     p.LineWidth = 0.95;
%     hold on
%     scatter(mean(Ea_store{n}), 0, 110, 'g', 'filled')
%     str = join(['Ea_', num2str(n)]);
%     title(str, 'FontSize',17)
%     if n==4
%         hold on
%         xline(24, 'Color', 'r', 'linewidth',3)
%         hold on
%         xline(36, 'Color', 'r', 'linewidth',3)
%     end
%     if n==1
%         hold on
%         xline(0, 'Color', 'r', 'linewidth',3)
%     end
% end
% 
% figure;
% for n = 1:4
%     p = subplot(2,2,n)
%     h = histogram(lnA_store{n});
%     h.FaceColor = [0,0.35,0.25];
%     h.EdgeColor = 'k'; %  [0.8, 0.8, 0.8];
%     h.LineWidth = 0.01;
%     h.FaceAlpha = 0.95;
%     p.LineWidth = 0.95;
%     hold on
%     scatter(mean(lnA_store{n}), 0, 110, 'g', 'filled')
%     str = join(['ln(A_', num2str(n),')']);
%     title(str, 'FontSize',17)
%     if n==4
%         hold on
%         xline(log(10^13.5), 'Color', 'r', 'linewidth',4)
%     end
% end




%% ARRHENIUS PLOTS
col{3} = [135, 171, 250]/256;
col{2} = [76, 117, 207]/256;
col{1} = [13, 48, 122]/256;
col{6} = [148, 38, 38]/256;
col{5} = [199, 48, 48]/256;
col{4} = [237, 164, 149]/256;

figure;
p =subplot(2,2,1);
p.LineWidth = 0.9;
scatter(1./T(idx), lnkk1(idx), 'k','filled')
hold on 
plot(1./T(idx), ln_kk1(idx),'Color',  col{3}, 'linewidth',2)
title('Adsorption k_1', 'FontSize', 15)
ylim([9.2,9.25])
set(gca, 'fontsize',13)
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

p =subplot(2,2,2)
p.LineWidth = 0.9;
scatter(1./T(idx), lnkk2(idx), 'k', 'filled')
hold on
plot(1./T(idx), ln_kk2(idx), 'Color', col{4}, 'LineWidth',2)
title('Adsorption k_2', 'FontSize', 15)
ylim([9.2,9.45])
set(gca, 'fontsize',13)
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

p =subplot(2,2,3)
p.LineWidth = 0.9;
scatter(1./T(idx), lnkk3(idx), 'k', 'filled')
hold on
plot(1./T(idx), ln_kk3(idx), 'Color', col{4}, 'LineWidth', 2)
title('Desorption k_3', 'FontSize', 15)
ylim([-6,1])
set(gca, 'fontsize',13)
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on

p = subplot(2,2,4)
p.LineWidth = 0.9;
scatter(1./T(idx), lnkk4(idx),  'k','filled')
hold on
plot(1./T(idx), ln_kk4(idx), 'Color', col{3}, 'LineWidth', 2)
title('Desorption k_4', 'FontSize', 15)
ylim([-10,0])
set(gca, 'fontsize',13)
xlabel('1/T [K^-^1]', 'FontSize', 15)
ylabel('ln(k)', 'FontSize', 15)
box on


%% ALL COVERAGE RESULTS
close all
col{3} = [135, 171, 250]/256;
col{2} = [76, 117, 207]/256;
col{1} = [13, 48, 122]/256;
col{6} = [148, 38, 38]/256;
col{5} = [199, 48, 48]/256;
col{4} = [237, 164, 149]/256;
% figure;
% for n = 1:N
%     blue = blue - [-0.07, 0.08, 0.09];
%     plot(time_area{n}(1:length(theta{n})),theta{n}, 'linewidth',1.8, 'Color', col{n})
%     hold on
% end
% set(gca, 'fontsize',15)
% xlabel('Time [s]','FontSize', 20)
% ylabel('Covereage [ML]', 'FontSize', 20)
% legend(temps_strings{1:end}, 'FontSize', 18)
% title('Inferred Latent States', 'FontSize', 20)
% grid on

% 
% 
%% MARKOV CHAINS 

% col{6} = [207, 37, 37]/256;
% 
% figure;
% p =subplot(2,2,1)
% p.LineWidth = 0.9;
% plot(x14chain(1:I,1), 'Color',col{1}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Adsorption k_1', 'FontSize', 15)
% xlim([0,I])
% box on
% 
% p = subplot(2,2,4)
% p.LineWidth = 0.9;
% plot(x14chain(1:I,2), 'Color',col{1},'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Desorption k_4', 'FontSize', 15)
% ylim([0,0.06])
% xlim([0,I])
% box on
% 
% p =subplot(2,2,2)
% p.LineWidth = 0.9;
% plot(x23chain(1:I,1), 'Color', col{6}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Adsorption k_2 ', 'FontSize', 15)
% xlim([0,I])
% box on
% 
% p = subplot(2,2,3)
% p.LineWidth = 0.9;
% plot(x23chain(1:I,2), 'Color', col{6}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Desorption k_3 ', 'FontSize', 15)
% ylim([0,1])
% xlim([0,I])
% box on
% 

%% ALL MCMC chains

% col{6} = [207, 37, 37]/256;
% 
% t = 6;
% 
% figure;
% p =subplot(2,2,1)
% p.LineWidth = 0.9;
% plot(kk1(t, 1:I), 'Color',col{1}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Adsorption k_1', 'FontSize', 15)
% xlim([0,I])
% box on
% 
% p = subplot(2,2,4)
% p.LineWidth = 0.9;
% plot(kk4(t, 1:I), 'Color',col{1},'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Desorption k_4', 'FontSize', 15)
% ylim([0,0.06])
% xlim([0,I])
% box on
% 
% p =subplot(2,2,2)
% p.LineWidth = 0.9;
% plot(kk2(t,1:I), 'Color', col{6}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Adsorption k_2 ', 'FontSize', 15)
% xlim([0,I])
% box on
% 
% p = subplot(2,2,3)
% p.LineWidth = 0.9;
% plot(kk3(t,1:I), 'Color', col{6}, 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Desorption k_3 ', 'FontSize', 15)
% ylim([0,1])
% xlim([0,I])
% box on



%% MIX figure for paper
% 
% figure;
% subplot(2,2,1)
% plot(k2chain(1:I), 'k', 'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Adsorption k_2 ', 'FontSize', 15)
% ylabel('Samples', 'FontSize', 15)
% xlabel('Iterations', 'FontSize', 15)
% xlim([0,I])
% 
% subplot(2,2,2)
% plot(k4chain(1:I), 'Color',col{1},'linewidth', 1)
% set(gca, 'fontsize',13)
% title('Desorption k_4', 'FontSize', 15)
% ylabel('Samples', 'FontSize', 15)
% xlabel('Iterations', 'FontSize', 15)
% ylim([0,0.03])
% xlim([0,I])
% 
% subplot(2,2,3)
% scatter(1./T(idx), lnkk2(idx), 'k', 'filled')
% hold on
% plot(1./T(idx), ln_kk2(idx), 'Color', col{4}, 'LineWidth',2)
% title('Adsorption k_2', 'FontSize', 15)
% ylim([9.1,9.6])
% set(gca, 'fontsize',13)
% xlabel('1/T [K^-^1]', 'FontSize', 15)
% ylabel('ln(k)', 'FontSize', 15)
% box on
% 
% subplot(2,2,4)
% scatter(1./T(idx), lnkk4(idx),  'k','filled')
% hold on
% plot(1./T(idx), ln_kk4(idx), 'Color', col{1}, 'LineWidth', 2)
% title('Desorption k_4', 'FontSize', 15)
% ylim([-10,-1])
% set(gca, 'fontsize',13)
% xlabel('1/T [K^-^1]', 'FontSize', 15)
% ylabel('ln(k)', 'FontSize', 15)
% box on

%% Variance convergences

figure;
for n = 1:N
 
    plot(I0:I, sig_store{n}(I0:end), 'linewidth',1.8, 'Color', col{n})
    hold on
end
set(gca, 'fontsize',15)
xlabel('Time [s]','FontSize', 20)
ylabel('Covereage [ML]', 'FontSize', 20)
legend(temps_strings{1:end}, 'FontSize', 18)
title('Inferred Latent States', 'FontSize', 20)
grid on





%% SAVE RESULTS
%save('Ea_res.mat', 'T', 'lnkk', 'ln_kk', "kk", 'R', 'Ea', 'ln_A', 'Ea_mean', "lnA_mean" )
