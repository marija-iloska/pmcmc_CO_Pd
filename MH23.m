function [x] = MH23(y, regions, x_old, bounds, P, dt, alpha3, alpha, alpha2, beta2)

% Restore last values
[~, ~, cov_sat, M, ~] = bounds{:};

% Old values
k2_old = x_old(1);
k3_old = x_old(2);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta3_old = (1 - k3_old*dt)*alpha3/(dt*k3_old);
 
% Propose sample for k4
k3_star = betarnd(alpha3, beta3_old)/dt;
a3_star = 1 - k3_star*dt;
a3_old = 1 - k3_old*dt;


% Propose sample for k1
x2_star = gamrnd(alpha2, beta2);
k2_star = k3_star*cov_sat/(M - cov_sat)/P + x2_star;
k2_lim = (M - a3_star*cov_sat)/(dt*P*(M - cov_sat));
k2_star = min([k2_lim, k2_star]);


% Get new porposal parameters
beta3_star = (1 - k3_star*dt)*alpha3/(k3_star*dt);


% Data
yR2 = y(regions{2});
yR2A = M - yR2;
yR3 = y(regions{3});


% Form a params
a2_star = k2_star*P*dt;
a2_old = k2_old*P*dt;

% Package
x_star = [k2_star, k3_star, x2_star];


% REGION III --------------------------------------------------------

% Mean of likelihood
mu_star =  2*a3_star*yR3(1:end-1);
mu_old = 2*a3_old*yR3(1:end-1);

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
ln_q3 = (beta3_star - 1)*log(1 - k3_old*dt) - (beta3_old - 1)*log(1 - k3_star*dt);

% Log likelihood ratio
ln_l3 = sum( (beta_star - beta_old).*log(1 - 2*yR3(2:end)));



% REGION II --------------------------------------------------------

% Mean of likelihood
mu_star =  2*(a2_star*yR2A(1:end-1) + a3_star*yR2(1:end-1));
mu_old = 2*(a2_old*yR2A(1:end-1) + a3_old*yR2(1:end-1));

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Log likelihood ratio
ln_l2 = sum((beta_star - beta_old).*(log(1 - 2*yR2(2:end))));



% ACCEPTANCE RATIO --------------------------------------------------
AR = exp(ln_l2 + ln_l3 + ln_q3); % + ln_q2);

% Accept or Reject
if (rand < AR)
    x = x_star;
else
    x = x_old;
end


end