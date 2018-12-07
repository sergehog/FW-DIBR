function [Sout] = recursive_compensation_weighted(S, W, sigma)

Sout = single(S);
Sout(isnan(S)) = 0;
Wx = single(~isnan(S)).*W;
Wx(isnan(Wx)) = 0;
S2 = recursive_gaussian(Sout.*Wx, sigma);
Wf = recursive_gaussian(Wx, sigma);

S2 = S2./Wf;
Sout(isnan(S)) = S2(isnan(S));