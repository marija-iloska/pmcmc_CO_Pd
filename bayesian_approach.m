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
    idx = J0:2:J;
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

% Ideal Gas constant  (kcal / (K mol))
R = 0.001987204258;

Regions = 4;

x = 1./(R*T);
N_samples = 2000;
mu = 10;
I0 = length(idx);

% Prior to the power I0-1
mu0 = mu/(I0-1);

% Local pos
for i = 1:I0

    y = {log(kk1(:,i)), log(kk2(:,i)), -log(kk3(:,i)), -log(kk4(:,i))};
    for r = 1:Regions
        temp = real(1/mu + sum((x'.*log(y{r}))));
        mu_local(r,i) = 1./temp;
    end
end

% Fused posterior
for r = 1:Regions
    temp = - 1/mu0 + sum(1./mu_local(r,:)/N^2);
    mu_f(r) = N/temp;
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

