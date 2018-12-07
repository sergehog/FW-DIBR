function [Iout, Zout] = rendering_hierarchical(Cv, I1, Z1, C1, radius, max_levels)
if nargin < 5
    radius = 3;
end
if nargin < 6
    [h, w] = size(Z1);
    max_levels = ceil(log2(min(h, w)));
end

[Iv, Zv, ~] = free_rendering(Cv, I1, Z1, C1);
if radius > 0
    Valid = imerode(~isnan(Zv),strel('disk',radius));
    Iv(repmat(~Valid, [1 1 3])) = nan;
    Zv(~Valid) = nan;
end
[Iout, Zout] = hierarchical_compensation_weighted(Iv, Zv, max_levels);
