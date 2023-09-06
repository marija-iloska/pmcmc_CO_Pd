function [x] = MH14op(y, regions, x_old, bounds, P, dt, var, alpha4, alpha, alpha1, beta1)

% Restore last values
[~, cut_off, M, ~] = bounds{:};

% Old values
k1_old = x_old(1);
k4_old = x_old(2);
x1_old = x_old(3);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta4_old = (1 - k4_old)*alpha4/k4_old;

% Propose sample for k4
k4_star = betarnd(alpha4, beta4_old);

var = 10;

% TRY LOG NORMAL
% mu_old_LN = 2;
% sigma2_old = 1;
% % sigma2_old = log( var*exp(-2/k4_old) + 1) + eps;
% % %sigma2_old = 2;
% % mu_old_LN = 1/k4_old - sigma2_old/2;
% lk4_star = lognrnd(mu_old_LN, sqrt(sigma2_old));
% k4_star = exp(- lk4_star)+eps;

% Compute parameters
% sigma2_star = log( var*exp(-2/k4_star) + 1) + eps;
% %sigma2_star = sigma2_old;
% mu_star_LN = 1/k4_star - sigma2_star/2;
% lk4_old = -log(k4_old);


% Propose sample for k1
% alpha1_old = x1_old^2/var;
% beta1_old = var/x1_old;
x1_star = gamrnd(alpha1, beta1);
k1_star = k4_star*cut_off/(M - cut_off)/P + x1_star;
k1_lim = (M - (1 - k4_star*dt)*cut_off)/(dt*P*(M - cut_off));

k1_star = min([k1_lim, k1_star]);

% Get new porposal parameters
beta4_star = (1 - k4_star)*alpha4/k4_star;
alpha1_star = x1_star^2/var;
beta1_star = var/x1_star;


% Data
yR1 = y(regions{1});
yR1A = M - yR1;
yR4 = y(regions{4});



% Form a params
a4_star = 1 - k4_star*dt;
a4_old = 1 - k4_old*dt;

a1_star = k1_star*P*dt;
a1_old = k1_old*P*dt;

% Package
x_star = [k1_star, k4_star, x1_star];



% REGION IV
mu_star =  2*a4_star*yR4(1:end-1);
mu_old = 2*a4_old*yR4(1:end-1);

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
% ln_q4 = log(lk4_star/lk4_old) + log(sigma2_old/sigma2_star) + 0.5*( (log(lk4_star) - mu_old_LN)^2/sigma2_old - ...
%     (log(lk4_old) - mu_star_LN)^2/sigma2_star );

ln_q4 = (beta4_star - 1)*log(1 - k4_old) - (beta4_old - 1)*log(1 - k4_star);

% Log likelihood ratio
ln_l4 = sum( (beta_star - beta_old).*log(1 - 2*yR4(2:end)));



% REGION I
mu_star =  2*(a1_star*yR1A(1:end-1) + a4_star*yR1(1:end-1));
mu_old = 2*(a1_old*yR1A(1:end-1) + a4_old*yR1(1:end-1));

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Log likelihood ratio
ln_l1 = sum((beta_star - beta_old).*(log(1 - 2*yR1(2:end))));

% Proposal ratios
% ln_q1 = log(gamma(x1_star)/gamma(x1_old)) + alpha1_old*log(beta1_old/x1_star) - alpha1_star*log(beta1_star/x1_old) +...
%     log(x1_star/x1_old) + x1_star/alpha1_old - x1_old/alpha1_star;
% 


% LOG Likelihood ratio
AR = exp(ln_l1 + ln_l4 + ln_q4); %+ ln_q1);

% Correction
if (isreal(AR) == 0)
    AR = 0;
end

% Accept or Reject
if (rand < AR)
    x = x_star;
else
    x = x_old;
end


end