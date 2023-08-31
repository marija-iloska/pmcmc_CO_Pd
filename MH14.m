function [x] = MH14(y, regions, a_old, bounds, alpha1, alpha4, beta1,  alpha)

% Restore last values
[~, cut_off, M, ~] = bounds{:};

% Old values
a1_old = a_old(1);
a4_old = a_old(2);

% Get old proposal for a4 q(a4 | alpha4, beta4_old)
beta4_old = (1 - a4_old)*alpha4/a4_old;


% Propose sample for a
a4_star = betarnd(alpha4, beta1);
y1_star = betarnd(alpha1, beta1);

%a4_star = min( abs(normrnd(0.5, 0.1)), 1);

% Get new proposal for a4 q(a4 | alpha4, beta4_star)
beta4_star = (1 - a4_star)*alpha4/a4_star;

% Data
yR1 = y(regions{1});
yR1A = M - yR1;
yR4 = y(regions{4});


% Max data point in R1
y1max = max(yR1);

% Upper bound of a1
a1_lim = max(0, (cut_off - a4_star*y1max)/(M - y1max));

% Sample a1
a1_star = a1_lim*y1_star;

% Concatenate params
a_star = [a1_star, a4_star];



% REGION IV
mu_star =  2*a4_star*yR4(1:end-1);
mu_old = 2*a4_old*yR4(1:end-1);

% Beta param
beta_star = (1-mu_star).*alpha./mu_star;
beta_old = (1-mu_old).*alpha./mu_old;

% Proposal ratios
ln_q4 = (beta4_star - 1)*log(1 - a4_old) - (beta4_old - 1)*log(1 - a4_star);

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



% LOG Likelihood ratio
AR = exp(ln_l1 + ln_l4); % + ln_q4);

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