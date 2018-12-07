function [Iv, Dv, Conf, Data] = rendering_dewith(Cv, I1, Z1, C1, fg_radius, patch, search, errode)

ZNear = min(Z1(:));
ZFar = max(Z1(:));
D1 = floor(255*(ZNear./Z1) .* ((ZFar-Z1)/(ZFar-ZNear)) + 0.5);
if 1==1
    fg_radius = 100;
    tic
    FgMask = local_otsu_thresholding(D1, fg_radius);
    toc
    %figure; imshow(FgMask, []); title('Local Otsu Thresholding'); drawnow;
else
    D1avg = fast_average(D1, fg_radius);
    FgMask = D1 > (D1avg+1);
    %clear D1avg D1;
end
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
[Iv, Dv, Conf, Data] = inpaint_rgbd_dewith(I0v, Depth, logical(FgMask), search, patch);

