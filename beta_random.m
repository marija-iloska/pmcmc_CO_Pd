function [Y] = beta_random(alpha, mu)

beta  = (1-mu).*alpha./mu;

Y = betarnd(alpha, beta);


end