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
    %str =join(['RESULTS/pmcmc_', temps_strings{n},'K_I3000', '.mat']);
    load(str)
    I0 = 2000;
    idx = I0:1:J;
    % Store chains
    kk1(n,:) = x14chain(idx,1);
    kk2(n,:) = x23chain(idx,1);
    kk3(n,:) = x23chain(idx,2);
    kk4(n,:) = x14chain(idx,2);

%     % Store Estimates
%     k1_adsorb(n) = k1_est;
%     k2_adsorb(n) = k2_est;
% 
%     k3_desorb(n) = k3_est;
%     k4_desorb(n) = k4_est;
    

end

% k1 = k1_adsorb;
% k2 = k2_adsorb;
% k3 = k3_desorb;
% k4 = k4_desorb;


clearvars -except kk1 kk2 kk3 kk4 k1 k2 k3 k4 J0 J idx
load Data/temps_info.mat
load Data/colors.mat

% Temperatures
T = [450, 460, 470, 475, 480, 490];

tidx = setdiff(1:N, []);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;
I0 = length(idx);

 %% NEW GAMMA
alpha = [0.1,0.1, 1, 1];
beta0 = [10,10, 10,10]*2;
ln_beta0 = [1, 2, 1, 1]*2;
ln_alpha = [1, 1, 10, 10];
 
% alpha = [0.0345,0.7242, 22.328, 30.0046];
% beta0 = [1, 1, 1, 1];
% ln_beta0 = [1, 1, 1, 1];
% ln_alpha = [9.2557, 10.0897, 21.3154, 24.4237];
x = 1./(R*T);
N_samples = 1;


for n = 1:4

    for i = 1:I0
    
        y = {log(kk1(:,i)), log(kk2(:,i)), -log(min(kk3(:,i),0.05)), -log(min(kk4(:,i),0.99))};
        beta_post = 1/( 1./beta0(n) + 1./sum(x'.*log(y{n})) ); 
        Ea_sample = gamrnd(alpha(n), beta_post, 1, N_samples);
        if (isreal(Ea_sample)==0)
            disp('stop')
        end
        Ea_store(i) = mean(Ea_sample);
    
        ln_beta_post = 1/(1/ln_beta0(n) - 1./sum(log(y{n})));
        lnA_sample = gamrnd(ln_alpha(n), ln_beta_post, 1, N_samples);
        lnA_store(i) = mean(lnA_sample);
    
    
    end
    
    nan_idx = find(isnan(Ea_store)==1);
    Ea_store(nan_idx) = [];
    Ea{n} = Ea_store;
    lnA{n} = lnA_store;

end

Ea_bayes = Ea;
lnA_bayes = lnA;

%save('RESULTS/bayes.mat', 'Ea_bayes', "lnA_bayes");


%%PLOT
dg = [29, 125, 80]/256;
figure;
for n = 1:4
    p = subplot(2,2,n);
    h = histogram(Ea{n});
    %h.FaceColor = [0,0,0.75];
    h.FaceColor = [186, 218, 247]/256;
    h.EdgeColor = [0.8, 0.8, 0.8];
    h.FaceAlpha = 0.67;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(Ea{n}), 0, 110, dg, 'filled')
    str = join(['Ea_', num2str(n)]);
    title(str, 'FontSize',17)
    if n==4
        hold on
        xline(24, 'Color', 'r', 'linewidth',3)
        hold on
        xline(36, 'Color', 'r', 'linewidth',3)
    end
    if n==1
        hold on
        xline(0, 'Color', 'r', 'linewidth',3)
    end
    %sgtitle('Gamma Prior')
    mean(Ea{n})
    box on
end

figure;
for n = 1:4
    p = subplot(2,2,n)
    h = histogram(lnA{n});
    h.FaceColor = [0.72,0.82,0.72];
    h.EdgeColor =  [0.8, 0.8, 0.8];
    h.FaceAlpha = 0.7;
    p.LineWidth = 0.95;
    hold on
    scatter(mean(lnA{n}), 0, 110, dg, 'filled')
    str = join(['ln(A_', num2str(n),')']);
    title(str, 'FontSize',17)
    if n==4
        hold on
        xline(log(10^13.5), 'Color', 'r', 'linewidth',3)
    end
    %sgtitle('Gamma Prior')
    mean(lnA{n})
end



