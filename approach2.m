close all
clear all
clc

load data/temps_info.mat


% Extract sampled parameters for each temperature
for n = 1 : N

    % Create filename and load data
    str = join(['results/', temps_strings{n},'K_J7000', '.mat']);

    load(str)

    I0 = 2000;
    idx = I0:1:J;
    
    % Store chains
    k1_chain(n,:) = x14chain(idx,1);
    k2_chain(n,:) = x23chain(idx,1);
    k3_chain(n,:) = x23chain(idx,2);
    k4_chain(n,:) = x14chain(idx,2);
  

end

clearvars -except k1_chain k2_chain k3_chain k4_chain J0 J idx N
load data/colors.mat

% Temperatures
T = [450, 460, 470, 475, 480, 490];

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;
I0 = length(idx);

 %% NEW GAMMA
% CASE 1
% alpha = [0.1,0.1, 1, 1];
% beta0 = [10,10, 10,10]*2;
% ln_beta0 = [1, 2, 1, 1]*2;
% ln_alpha = [1, 1, 10, 10];

% CASE 2 and CASE 3 (adjust)
alpha = [0.0345,0.7242, 22.328, 30.0046]/1;
beta0 = [1, 1, 1, 1]*1;
ln_beta0 = [1, 1, 1, 1]*3;
ln_alpha = [9.2557, 10.0897, 21.3154, 24.4237]/3;
x = 1./(R*T);
N_samples = 1;

% Loop for each region
for r = 1:4

    % Compute posterior for every MCMC sample
    for i = 1:I0
    
        % Transform data
        y = {log(k1_chain(:,i)), log(k2_chain(:,i)), -log(min(k3_chain(:,i),0.05)), -log(min(k4_chain(:,i),0.99))};

        % Compute Ea posterior
        beta_post = 1/( 1./beta0(r) + 1./sum(x'.*log(y{r})) ); 
        Ea_sample = gamrnd(alpha(r), beta_post, 1, N_samples);
        Ea_store(i) = mean(Ea_sample);

        % Compute lnA posterior
        ln_beta_post = 1/(1/ln_beta0(r) - 1./sum(log(y{r})));
        lnA_sample = gamrnd(ln_alpha(r), ln_beta_post, 1, N_samples);
        lnA_store(i) = mean(lnA_sample);
        
    end
    
    % Store means
    Ea_approach2{r} = Ea_store;
    lnA_approach2{r} = lnA_store;

end

%save('results/approach2.mat', 'Ea_approach2', "lnA_approach2");


%% VISUALIZE RESULTS 

% Choose color
dg = [29, 125, 80]/256;

% Create histograms for Ea
figure;
for r = 1:4
    p = subplot(2,2,r);
    h = histogram(Ea_approach2{r});
    %h.FaceColor = [0,0,0.75];
    h.FaceColor = [186, 218, 247]/256;
    h.EdgeColor = [0.8, 0.8, 0.8];
    h.FaceAlpha = 0.67;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(Ea_approach2{r}), 0, 110, dg, 'filled')
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
    box on
end

% Create histograms for lnA
figure;
for r = 1:4
    p = subplot(2,2,r);
    h = histogram(lnA_approach2{r});
    h.FaceColor = [0.72,0.82,0.72];
    h.EdgeColor =  [0.8, 0.8, 0.8];
    h.FaceAlpha = 0.7;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(lnA_approach2{r}), 0, 110, dg, 'filled')
    str = join(['ln(A_', num2str(r),')']);
    title(str, 'FontSize',17)
    if r==4
        hold on
        xline(log(10^13.5), 'Color', 'r', 'linewidth',3)
    end
end

