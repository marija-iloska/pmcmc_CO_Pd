function [Ea, ln_A, ln_k_fit] = compute_Ea(lnk, T, T_fit, R)


% Activation energy - Ea
% The units of Ea depend on the units of the R value entered

% Arrhenious equation
% k = A exp(-Ea/(RT))
% ln(k) = -Ea/(RT) + ln(A)

% Let x = - 1/(RT)
% ln(k) = Ea x + ln(A)
x = - 1./(R*T);
x_fit = -1./(R*T_fit);

% Linear Model fitting
dlm_Ea = fitlm(x, lnk,'Intercept',true);
Ea = dlm_Ea.Coefficients.Estimate(2);
ln_A = dlm_Ea.Coefficients.Estimate(1);

% Get fitted ln(k)s
ln_k_fit = Ea*x_fit + ln_A;

end