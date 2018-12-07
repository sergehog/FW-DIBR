function [Ifg, Zfg, Ibg, Zbg, Alpha] = recursive_decomposition(I, Zv, avg_alpha, erode)
if nargin < 4
    erode = 1;
end
[c] = size(I, 3);
ZNear = min(Zv(:));
ZFar = max(Zv(:));
Depth = floor(255*(ZNear./Zv) .* ((ZFar-Zv)/(ZFar-ZNear)) + 0.5);
W = single(~isnan(Depth));
Depth(isnan(Depth)) = 0;
Davg = recursive_gaussian(Depth, avg_alpha)./recursive_gaussian(W, avg_alpha);

FgMask = Depth > (Davg + 2) & (W>0);
%FgMask2 = Depth > (Davg + 3) & (W>0);
%CMask = (FgMask > FgMask2) & (W>0);
BgMask = (~(Depth > (Davg + 4))) & (W > 0);
BgMask = imerode(BgMask, strel('disk',erode)); 

%[Dx, Dy] = gradient(D);
%Edges = sqrt(Dx.^2 + Dy.^2) > 0.5;
%Edges = imdilate(Edges, strel('disk',2)) & imdilate(BgMask, strel('disk',2));
%BgMask = int8(BgMask) - int8(Edges) > 0;


Ifg = I;
Ifg(repmat(~FgMask, [1 1 c])) = nan;

Zfg = Zv;
Zfg(~FgMask) = nan;

Ibg = I;
Ibg(repmat(~BgMask, [1 1 c])) = nan;

Zbg = Zv;
Zbg(~BgMask) = nan;

Alpha = single(FgMask);

