function [Iout, Zout] = recursive_compensation(Iv, Zv, sigma)
maxZ = max(Zv(:));
minZ = min(Zv(:));

Iout = single(Iv);
Iout(isnan(Iv)) = 0;
%W = single(~isnan(Zv));
W = (Zv-minZ)/(maxZ-minZ);
W(isnan(W)) = 0;
Iout = Iout.*repmat(W, [1 1 3]);
Zout = Zv.*W;
Zout(isnan(Zout)) = 0;

If = recursive_gaussian(Iout, sigma);
Zf = recursive_gaussian(Zout, sigma);
Wf = recursive_gaussian(W, sigma);

If = If./repmat(Wf,[1 1 3]);
Zf = Zf./Wf;
Iout = single(Iv);
Iout(isnan(Iv)) = If(isnan(Iv));
Zout = single(Zv);
Zout(isnan(Zv)) = Zf(isnan(Zv));