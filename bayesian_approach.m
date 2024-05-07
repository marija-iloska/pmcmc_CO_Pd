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
    J0 = 2000;
    idx = J0:4:J;
    % Store chains
    kk1(n,:) = x14chain(idx,1);
    kk2(n,:) = x23chain(idx,1);
    kk3(n,:) = x23chain(idx,2);
    kk4(n,:) = x14chain(idx,2);

end

clearvars -except kk1 kk2 kk3 kk4 J0 J idx
load Data/temps_info.mat
load Data/colors.mat

% Temperatures
T = [450, 460, 470, 475, 480, 490];

tidx = setdiff(1:N, []);

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;

Regions = 4;

x = 1./(R*T(tidx));
N_samples = 2000;
mu = [100, 100, 0.001, 0.001];
I0 = length(idx);



% Set negative for R12
y = log(kk4);

% Priors 
lambda_Ea = 50;
lambda_lnA = 10;





for i = 1:I0

    % Sample from prior
    M = 30;
    Ea = exprnd(lambda_Ea, 1,M);
    lnA = exprnd(lambda_lnA, 1,M);


    % Reusable likelihood terms
    term1 = sum(y(:,i),1);
    term1 = sum(term1);
    term2 = prod(y(:,i), 1);
    term3 = sum(y(:,i)./T')/R;


    term2 = lnA*log(term2);
    term3 = -Ea*log(-term3);

    logL = term1 + term3 + term2;


    % Weights
    log_w = logL;
    w = exp(log_w - max(log_w));

    % Final posterior
    Ea4_mean(i) = sum(w.*Ea);
    lnA4_mean(i) = sum(w.*lnA);

end



% Evaluate
w = fused_likelihood_weights(term1, term2, term3, Ea, lnA, lambda_Ea, lambda_lnA, I0);







% Prior to the power I0-1
mu0 = mu./(I0-1);

% Local pos
for i = 1:I0

    y = {log(kk1(tidx,i)), log(kk2(tidx,i)), -log(kk3(tidx,i)), -log(kk4(tidx,i))};
    for r = 1:Regions
%         temp = real(1/mu + sum((x'.*log(y{r}))));
%         mu_local(r,i) = 1./temp;
          lambda_local =  real(1./(log(R) - sum(log(y{r}./T(tidx)')) ));
          mu_local(r,i) = 1./lambda_local;
    end
end

% Fused posterior
for r = 1:Regions
    temp = sum(1./mu_local(r,:))/I0;
    mu_f(r) = mu0(r) + 1./temp;
end








for r = 1:Regions

    Ea{r} = exprnd(mu_f(r), 1, N_samples);
    subplot(2,2,r)
    hist(Ea{r})
    hold on
    scatter(mean(Ea{r}), 0, 70, 'g', 'filled')
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
    mean(Ea{r})
end
sgtitle('Exponential Prior')

