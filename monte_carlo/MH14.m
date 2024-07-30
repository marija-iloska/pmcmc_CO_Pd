function [x] = MH14(y, regions, x_old, bounds, P, dt, alpha4, alpha, alpha1, beta1)

% Restore last values
[~, cut_off, ~, M, ~] = bounds{:};

% Old values
k1_old = x_old(1);
k4_old = x_old(2);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta4_old = (1 - k4_old*dt)*alpha4/(k4_old*dt);

% Propose sample for k4
k4_star = betarnd(alpha4, beta4_old)/dt;
a4_star = 1 - k4_star*dt;
a4_old = 1 - k4_old*dt;

% Propose sample for k1
x1_star = gamrnd(alpha1, beta1);
k1_star = k4_star*cut_off/(M - cut_off)/P + x1_star;
k1_lim = (M - a4_star*cut_off)/(dt*P*(M - cut_off));
k1_star = min([k1_lim, k1_star]);

% Get new porposal parameters
beta4_star = (1 - k4_star*dt)*alpha4/(k4_star*dt);



% Data
yR1 = y(regions{1});
yR1A = M - yR1;
yR4 = y(regions{4});



% Form a params
a1_star = k1_star*P*dt;
a1_old = k1_old*P*dt;

% Package
x_star = [k1_star, k4_star, x1_star];


% REGION IV --------------------------------------------------------

% Mean of likelihood
mu_star =  2*a4_star*yR4(1:end-1);
mu_old = 2*a4_old*yR4(1:end-1);

% Beta param old and proposed 
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
ln_q4 = (beta4_star - 1)*log(1 - k4_old*dt) - (beta4_old - 1)*log(1 - k4_star*dt);

% Log likelihood ratio
ln_l4 = sum( (beta_star - beta_old).*log(1 - 2*yR4(2:end)));



% REGION I--------------------------------------------------------

% Mean of likelihood
mu_star =  2*(a1_star*yR1A(1:end-1) + a4_star*yR1(1:end-1));
mu_old = 2*(a1_old*yR1A(1:end-1) + a4_old*yR1(1:end-1));

% Beta param old and proposed
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Log likelihood ratio
ln_l1 = sum((beta_star - beta_old).*(log(1 - 2*yR1(2:end))));


% ACCEPTANCE RATIO --------------------------------------------------
AR = exp(ln_l1 + ln_l4 + ln_q4);

% Accept or Reject
if (rand < AR)
    x = x_star;
else
    x = x_old;
end


end