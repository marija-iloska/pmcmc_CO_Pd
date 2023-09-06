function [x] = MH23op(y, regions, x_old, bounds, P, dt, var, alpha3, alpha, alpha2, beta2, cov_sat)

% Restore last values
[~, cut_off, M, ~] = bounds{:};




% Old values
k2_old = x_old(1);
k3_old = x_old(2);
x2_old = x_old(3);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta3_old = (1 - k3_old)*alpha3/k3_old;
 
 % Propose sample for k4
k3_star = betarnd(alpha3, beta3_old);

var = 10;

% TRY LOG NORMAL
% mu_old_LN = 2;
% sigma2_old = 1;
% sigma2_old = log( var*exp(-2/k4_old) + 1) + eps;
% %sigma2_old = 2;
% mu_old_LN = 1/k4_old - sigma2_old/2;
% lk3_star = lognrnd(mu_old_LN, sqrt(sigma2_old));
% k3_star = exp(- lk3_star)+eps;



% Propose sample for k1
% alpha2_old = x2_old^2/var;
% beta2_old = var/x2_old;
x2_star = gamrnd(alpha2, beta2);
k2_star = k3_star*cov_sat/(M - cov_sat)/P + x2_star;
k2_lim = (M - (1 - k3_star*dt)*cov_sat)/(dt*P*(M - cov_sat));

k2_star = min([k2_lim, k2_star]);


% Get new porposal parameters
beta3_star = (1 - k3_star)*alpha3/k3_star;
alpha2_star = x2_star^2/var;
beta2_star = var/x2_star;


% Data
yR2 = y(regions{2});
yR2A = M - yR2;
yR3 = y(regions{3});



% Form a params
a3_star = 1 - k3_star*dt;
a3_old = 1 - k3_old*dt;

a2_star = k2_star*P*dt;
a2_old = k2_old*P*dt;

% Package
x_star = [k2_star, k3_star, x2_star];


% REGION III
mu_star =  2*a3_star*yR3(1:end-1);
mu_old = 2*a3_old*yR3(1:end-1);

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
ln_q3 = (beta3_star - 1)*log(1 - k3_old) - (beta3_old - 1)*log(1 - k3_star);

% Log likelihood ratio
ln_l3 = sum( (beta_star - beta_old).*log(1 - 2*yR3(2:end)));



% REGION II
mu_star =  2*(a2_star*yR2A(1:end-1) + a3_star*yR2(1:end-1));
mu_old = 2*(a2_old*yR2A(1:end-1) + a3_old*yR2(1:end-1));

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Log likelihood ratio
ln_l2 = sum((beta_star - beta_old).*(log(1 - 2*yR2(2:end))));

% Proposal ratio
% ln_q2 = log(gamma(x2_star)/gamma(x2_old)) + alpha2_old*log(beta2_old/x2_star) - alpha2_star*log(beta2_star/x2_old) +...
%     log(x2_star/x2_old) + x2_star/alpha2_old - x2_old/alpha2_star;


% LOG Likelihood ratio
AR = exp(ln_l2 + ln_l3 + ln_q3); % + ln_q2);

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