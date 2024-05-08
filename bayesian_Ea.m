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



k1 = k1_adsorb;
k2 = k2_adsorb;
k3 = k3_desorb;
k4 = k4_desorb;





%% EXP PRIOR Ea
% clc
% x = 1./(R*T);
% N_samples = 10000;
% lambda = [1000, 1000, 0.1, 0.1];
% 
% 
% y = {log(k1), log(k2), -log(k3), -log(k4)};
% 
% figure;
% for n = 1:4
%     lambda_post = 1./( 1./lambda(n) + log(R) - sum(log(y{n}./T)) );
%     mu = 1/lambda_post;
%     Ea{n} = exprnd(mu, 1,N_samples);
%     subplot(2,2,n)
%     hist(Ea{n})
%     hold on
%     scatter(mean(Ea{n}), 0, 70, 'g', 'filled')
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
%     mean(Ea{n})
% end
% sgtitle('Exponential Prior')

%% OLD GAMMA prior Ea
% alpha = [2,2,0.1,0.1];
% beta_a = 100;
% x = -1./(R*T);
% 
% figure;
% for n = 1:4
%     beta_post = ((1 - beta_a*sum(x.*log(y{n})))/beta_a); 
%     Ea{n} = gamrnd(alpha(n), beta_post, 1, N_samples);
%     subplot(2,2,n)
%     hist(Ea{n})
%     hold on
%     scatter(mean(Ea{n}), 0, 70, 'g', 'filled')
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
% sgtitle('Gamma Prior')

%% PRE EXP FACTOR
% figure;
% for n = 1:4
%     lambda_post = 1/( 1/lambda + sum(log(y{n}))) ;
%     mu = 1/lambda_post;
%     lnA{n} = exprnd(mu, 1,N_samples);
%     subplot(2,2,n)
%     hist(lnA{n})
%     hold on
%     scatter(mean(lnA{n}), 0, 70, 'g', 'filled')
%     str = join(['ln(A_', num2str(n),')']);
%     title(str, 'FontSize',17)
%     if n == 4
%         hold on
%         xline(log(10^13.5), 'Color', 'r', 'linewidth',3)
%     end
% end
% sgtitle('Exponential Prior')
% 
% 

% %% NEW GAMMA
% alpha = 200;
% beta0 = 1;
% x = 1./(R*T);
% N_samples = 1000;
% 
% y = {log(k1), log(k2), -log(k3), -log(k4)};
% 
% figure;
% for n = 3
%     beta_post = 1/( 1/beta0 + sum(x.*log(y{n})) ); 
%     Ea{n} = gamrnd(alpha, beta_post, 1, N_samples);
%     subplot(2,2,n)
%     hist(Ea{n})
%     hold on
%     scatter(mean(Ea{n}), 0, 70, 'g', 'filled')
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
% sgtitle('Gamma Prior')


