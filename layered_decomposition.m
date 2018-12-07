function [Df, Db, Ifg, Ibg, alpha] = layered_decomposition(I1, D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg)

matting_radius = 1;
matting_sigma = 50;
matting_epsilon = 1e-4;
matting_lambda = 0.1;

depth_thr0 = 1;
depth_thr1 = 2;

[h, w] = size(D1);
HW = h * w;

D = D1;
D(isnan(D)) = 0;
D_hat = recursive_gaussian(D, sigma_depth)./recursive_gaussian(single(~isnan(D1)), sigma_depth);
clear D;
Mask0 = D1 > (D_hat + depth_thr0) & ~isnan(D1);
Mask1 = D1 > (D_hat + depth_thr1) & ~isnan(D1);
OverlapMask = Mask0 & ~Mask1;

FgDepthMask = imerode(Mask0 | OverlapMask, strel('disk', errosion_depth)) | OverlapMask;
BgDepthMask = imerode(~Mask1 | OverlapMask, strel('disk', errosion_depth)) | OverlapMask;

Df = D1;
Df(~FgDepthMask) = nan;
Db = D1;
Db(~BgDepthMask) = nan;


FgColorMask = imerode(Mask0 & ~isnan(D1), strel('disk', errosion_fg)) | OverlapMask;
%BgColorMask = (imerode((~Mask1 & ~isnan(D1)), strel('disk', errosion_bg)) | OverlapMask) & (~FgColorMask);
BgColorMask = imerode((~Mask1 & ~isnan(D1)), strel('disk', errosion_bg)) | OverlapMask;

Ti = FgColorMask + 1 - (BgColorMask & ~OverlapMask);
%figure; imshow(Ti, [0 2]); title('Trimap'); drawnow;

AL = get_learning_laplacian(uint8(I1), uint8(Ti), matting_radius, matting_sigma, matting_epsilon);
S = double(Ti(:)~=1);
S = spdiags(S(:), 0, HW, HW);
alpha = (matting_lambda*AL + S'*S) \ (S'*(double(Ti(:))));
alpha = single(reshape(alpha, [h w]))/2;
%alpha = single(reshape(alpha, [h w]) + 1)/2;
alpha(alpha < 0) = 0;
alpha(alpha > 1) = 1;
%Ifg = 0;
%Ibg = 0;
%return
%%

%[Ifg, Ibg] = solveFB(double(I1)/255, double(alpha));
[Ifg, ~] = solveFB(double(I1)/255, double(alpha));
Ifg = single(Ifg*255);
%Ibg = single(Ibg*255);    
Ifg(repmat(alpha<=0.1, [1 1 3])) = nan;
Ifg(repmat(OverlapMask, [1 1 3])) = I1(repmat(OverlapMask, [1 1 3]));    

%Ifg_f = Ifg;
%Ifg_f(isnan(Ifg)) = 0;
%Ifg_f = recursive_gaussian(Ifg_f, 0.8)./recursive_gaussian(single(~isnan(Ifg)), 0.8);
%Ifg(isnan(Ifg)) = Ifg_f(isnan(Ifg));
%clear Ifg_f;
%Ibg(repmat(MDiff, [1 1 3])) = I1(repmat(MDiff, [1 1 3]));

Ibg = I1;
Ibg(repmat(alpha > 0.02, [1 1 3])) = nan;
Ibg(repmat(~BgColorMask, [1 1 3])) = nan;
Ibg(repmat(OverlapMask, [1 1 3])) = I1(repmat(OverlapMask, [1 1 3]));    
%Ibg(repmat(~BgColorMask, [1 1 3])) = nan;
%Ibg(repmat(isnan(Db), [1 1 3])) = nan;
