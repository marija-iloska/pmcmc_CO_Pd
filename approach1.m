close all
clear all
clc

% Load data
load data/temps_info.mat
load data/colors.mat

% Extract sampled parameters for each temperature
for n = 1 : N

    % Create filename and load data
    str = join(['results/', temps_strings{n},'K_J7000', '.mat']);

    load(str)
    theta{n} = theta_est;

    % Store chains
    k1_chain(n,:) = x14chain(:,1);
    k2_chain(n,:) = x23chain(:,1);
    k3_chain(n,:) = x23chain(:,2);
    k4_chain(n,:) = x14chain(:,2);

end

% Clear some storage
clearvars -except k1_chain k2_chain k3_chain k4_chain J0 J idx N col time_area theta


%% GET ESTIMATES

% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Burn-in
I0 = 2000;

% In case we want to exclude a tempereature
idx = setdiff(1:N, []);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;

% Transform chain results to logs
lnk1_chain = log(k1_chain);
lnk2_chain = log(k2_chain);
lnk3_chain = log(k3_chain);
lnk4_chain = log(k4_chain);

T_fit = T(1):1:T(end);

% Compute Ea samples for each region using k samples (after burn-in)
for j = 1:(J-I0-1)
    [Ea4(j), ln_A4(j), lnk4_fit(:,j)] = compute_Ea(lnk4_chain(:,I0+j), T, T_fit, R);
    [Ea3(j), ln_A3(j),  lnk3_fit(:,j)] = compute_Ea(lnk3_chain(:,I0+j), T, T_fit, R);
    
    [Ea2(j), ln_A2(j),  lnk2_fit(:,j)] = compute_Ea(lnk2_chain(:,I0+j), T, T_fit, R);
    [Ea1(j), ln_A1(j),  lnk1_fit(:,j)] = compute_Ea(lnk1_chain(:,I0+j), T, T_fit, R);
end

% Compute mean estimates
Ea_mean = mean([Ea1; Ea2; Ea3; Ea4],2);
lnA_mean = mean([ln_A1; ln_A2; ln_A3; ln_A4],2);


% Store results
Ea_approach1 = {Ea1, Ea2, Ea3, Ea4};
lnA_approach1 = {ln_A1, ln_A2, ln_A3, ln_A4};


% Fitted lnk
lnk1_fit_est = mean(lnk1_fit, 2);
lnk2_fit_est = mean(lnk2_fit, 2);
lnk3_fit_est = mean(lnk3_fit, 2);
lnk4_fit_est = mean(lnk4_fit, 2);


% From data lnk
lnk1_est = mean(lnk1_chain,2);
lnk2_est = mean(lnk2_chain,2);
lnk3_est = mean(lnk3_chain,2);
lnk4_est = mean(lnk4_chain,2);

% For plotting
lnk_data = {lnk1_est, lnk2_est, lnk3_est, lnk4_est};
lnk_fit = {lnk1_fit_est, lnk2_fit_est, lnk3_fit_est, lnk4_fit_est};

% Clear some storage
clearvars lnk1_chain lnk2_chain lnk3_chain lnk4_chain 

% SAVE RESULTS
%save('results/approach1.mat', 'Ea_approach1', "lnA_approach1");


%% PARAMETER HISTOGRAMS 

% ACTIVATION ENERGY  --------------------------------------------
figure;
for r = 1:4
    p = subplot(2,2,r);
    h = histogram(Ea_approach1{r});
    h.FaceColor = [0, 0.35, 0.65];
    h.EdgeColor = 'k';
    h.LineWidth = 0.01;
    h.FaceAlpha = 0.95;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(Ea_approach1{r}), 0, 110, 'g', 'filled')
    str = join(['Ea_', num2str(r)]);
    title(str, 'FontSize',17)
    if r==4
        hold on
        xline(24, 'Color', 'r', 'linewidth',3)
        hold on
        xline(36, 'Color', 'r', 'linewidth',3)
    end
    if r==1
        hold on
        xline(0, 'Color', 'r', 'linewidth',3)
    end
end

% PRE-EXPONENTIAL FACTORS  --------------------------------------------
figure;
for r = 1:4
    p = subplot(2,2,r);
    h = histogram(lnA_approach1{r});
    h.FaceColor = [0,0.35,0.25];
    h.EdgeColor = 'k'; 
    h.LineWidth = 0.01;
    h.FaceAlpha = 0.95;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(lnA_approach1{r}), 0, 110, 'g', 'filled')
    str = join(['ln(A_', num2str(r),')']);
    title(str, 'FontSize',17)
    if r==4
        hold on
        xline(log(10^13.5), 'Color', 'r', 'linewidth',4)
    end
end


%% ARRHENIUS PLOTS  ---------------------------------------------

figure;
for r = 1:4
    p =subplot(2,2,r);
    p.LineWidth = 0.9;
    scatter(1./T(idx), lnk_data{r}, 'k','filled')
    hold on 
    plot(1./T_fit, lnk_fit{r},'Color',  col{r}, 'linewidth',2)
    title('Adsorption k_1', 'FontSize', 15)
    set(gca, 'fontsize',13)
    xlabel('1/T [K^-^1]', 'FontSize', 15)
    ylabel('ln(k)', 'FontSize', 15)
    box on
end



% ALL COVERAGE RESULTS  -------------------------------------------------------
blue = [1,1,1];
figure;
for n = 1:N
    blue = blue - [-0.07, 0.08, 0.09];
    plot(time_area{n}(1:length(theta{n})),theta{n}, 'linewidth',1.8, 'Color', col{n})
    hold on
end
set(gca, 'fontsize',15)
xlabel('Time [s]','FontSize', 20)
ylabel('Covereage [ML]', 'FontSize', 20)
legend(temps_strings{1:end}, 'FontSize', 18)
title('Inferred Latent States', 'FontSize', 20)
grid on


% ALL MCMC chains  ---------------------------------------------------------------------

k_chains = {k1_chain, k2_chain, k3_chain, k4_chain};
chain_plot_labels = {'Adsorption k_1', 'Adsorption k_2', 'Desorption k_3', 'Desorption k_4'};

% Choose temperature
n = 6;

figure;
for r = 1:4
    subplot(2,2,r)
    plot(k_chains{r}(n,1:J), 'Color',col{r}, 'linewidth', 1)
    title(chain_plot_labels{r}, 'FontSize', 15)
    xlim([0,J])
end

