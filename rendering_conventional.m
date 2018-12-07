function [Iout, Zv] = rendering_conventional(Cv, I1, Z1, C1, alpha_hf, radius)
[Iv, Zv, ~] = free_rendering(Cv, I1, Z1, C1);    
if radius > 0
    %Valid = imopen(~isnan(Zv),strel('disk',radius));
    Valid = imerode(~isnan(Zv),strel('disk',radius));
    Iv(repmat(~Valid, [1 1 3])) = nan;
    Zv(~Valid) = nan;
end
[Iout] = recursive_compensation(Iv, Zv, alpha_hf);
