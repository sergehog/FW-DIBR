function [Iv, Dv, Conf, Data] = rendering_olsson(Cv, I1, Z1, C1, fg_radius, patch, search, errode)
[h, w, ~] = size(I1);
[UV, V] = meshgrid(0:(w-1), 0:(h-1));
UV(:,:,2) = V;
clear V;
UV = single(UV);

ZNear = min(Z1(:));
ZFar = max(Z1(:));
D1 = floor(255*(ZNear./Z1) .* ((ZFar-Z1)/(ZFar-ZNear)) + 0.5);

D1avg = fast_average(D1, fg_radius);
FgMask = D1 > (D1avg+1);
clear D1avg D1;

[I0v, Z0v, ~] = free_rendering(Cv, I1, Z1, C1);
[FgMask, ~, ~] = free_rendering(Cv, single(FgMask), Z1, C1);
Valid = imerode(~isnan(Z0v), strel('disk', errode));
I0v(repmat(~Valid, [1 1 3])) = nan;
Depth = floor(255*(ZNear./Z0v) .* ((ZFar-Z0v)/(ZFar-ZNear)) + 0.5);
Depth(~Valid) = nan;
FgMask(isnan(FgMask)) = 0;
%figure; imshow(uint8(I0v)); title('Rendered Image with Holes');
%figure; imshow(Depth, [0 255]); colormap(pink); title('Rendered Depth with Holes');
%figure; imshow(FgMask, [0 1]); title('FgMask');
%drawnow;
[Iv, Dv, Conf, Data] = inpaint_rgbd_olsson(I0v, Depth, FgMask>0.1, search, patch);

