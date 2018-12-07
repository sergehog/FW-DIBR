function [Iv, Dv, Conf, Data] = rendering_jtdi(Cv, I1, Z1, C1, search, inpaint, errode)
[I0v, Z0v, ~] = free_rendering(Cv, I1, Z1, C1);
ZNear = min(Z1(:));
ZFar = max(Z1(:));
Depth = floor(255*(ZNear./Z0v) .* ((ZFar-Z0v)/(ZFar-ZNear)) + 0.5);
Valid = imerode(~isnan(Z0v), strel('disk', errode));
I0v(repmat(~Valid, [1 1 3])) = nan;
Depth(~Valid) = nan;
%figure; imshow(uint8(I0v)); title('Rendered Image with Holes');
%figure; imshow(Depth, [0 255]); colormap(pink); title('Rendered Depth with Holes');
%drawnow;
[Iv, Dv, Conf, Data] = inpaint_rgbd_jtdi(I0v, Depth, search, inpaint);
