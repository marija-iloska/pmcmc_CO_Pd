function [Ea, A, Ea_SE, A_SE, ln_k, Rsq_Ea] = get_Ea(k, T, R)

% The units of Ea depend on the units of the R value entered

% Arrhenious equation
% k = A exp(-Ea/(RT))
% ln(k) = -Ea/(RT) + ln(A)

% Let x = - 1/(RT)
% ln(k) = Ea x + ln(A)
x = - 1./(R*T);

% Linear Model fitting
dlm_Ea = fitlm(x, log(k),'Intercept',true);
Ea = dlm_Ea.Coefficients.Estimate(2);
Ea_SE = dlm_Ea.Coefficients.SE(2);
ln_A = dlm_Ea.Coefficients.Estimate(1);
A = exp(ln_A);
A_SE = dlm_Ea.Coefficients.SE(1);

% R squared

% Get fitted ln(k)s
ln_k = Ea*x + ln_A;
Rsq_Ea = dlm_Ea.Rsquared.Ordinary;


end