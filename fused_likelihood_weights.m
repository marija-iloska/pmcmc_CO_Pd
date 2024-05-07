function [w] = fused_likelihood_weights(term1, term2, term3, Ea, lnA, lambda_Ea, lambda_lnA, I0)


term2 = lnA*sum(log(term2));
term3 = -Ea*sum(log(-term3));

% Loglikelihood
logL = -term1 - term3 - term2;

% Priors
logq = -Ea/lambda_Ea - lnA/lambda_lnA;

% Weights
log_w = logL - logq;
w = exp(log_w - max(log_w));

end