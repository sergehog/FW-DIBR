function [Sout] = recursive_compensation_simple(S, sigma)

Sout = single(S);
Sout(isnan(S)) = 0;
W = single(~isnan(S));

S2 = recursive_gaussian(Sout, sigma);
Wf = recursive_gaussian(W, sigma);

S2 = S2./Wf;
Sout(isnan(S)) = S2(isnan(S));