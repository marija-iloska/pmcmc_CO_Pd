function [x] = MH23(y, regions, a_old, bounds, alpha2, alpha3, beta2,  alpha, cov_sat)

% Restore last values
[~, cut_off, M, ~] = bounds{:};

% Old values
a2_old = a_old(1);
a3_old = a_old(2);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta3_old = (1 - a3_old)*alpha3/a3_old;


% Propose sample for a
a3_star = betarnd(alpha3, beta2);
y2_star = betarnd(alpha2, beta2);

%a3_star = min( abs(normrnd(0.5, 0.1)), 1);

% Get new proposal for a4 q(a4 | alpha4, beta4_star)
beta3_star = (1 - a3_star)*alpha3/a3_star;

% Data
yR2 = y(regions{2});
yR2A = M - yR2;
yR3 = y(regions{3});


% Max data point in R1
y2max = max(yR2);

% Upper bound of a1
a2_lim = (M - a3_star*cov_sat)/(M - cov_sat);

% Sample a1
a2_star = a2_lim*y2_star;

% Concatenate params
a_star = [a2_star, a3_star];



% REGION III
mu_star =  2*a3_star*yR3(1:end-1);
mu_old = 2*a3_old*yR3(1:end-1);

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
ln_q3 = (beta3_star - 1)*log(1 - a3_old) - (beta3_old - 1)*log(1 - a3_star);

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



% LOG Likelihood ratio
AR = exp(ln_l2 + ln_l3); % + ln_q3);

% Correction
if (isreal(AR) == 0)
    AR = 0;
end

% Accept or Reject
if (rand < AR)
    x = a_star;
else
    x = a_old;
end


end