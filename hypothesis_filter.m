function D_hat = hypothesis_filter(D, I, sigma_color, alpha_spacial, cost_thr)
if nargin < 3
    sigma_color = 0.01;    
end
if nargin < 4
    alpha_spacial = 0.88;    
end
if nargin < 5
    cost_thr = 3;    
end

maxdisp = ceil(max(D(:)));
mindisp = floor(min(D(:)));
layers = maxdisp - mindisp + 1;
[h, w] = size(D);

Cost = zeros([h w layers+1], 'single');
for i=1:layers
    d = i + mindisp - 1;
    Cost(:,:,i) = single(min(cost_thr, abs(D-d)));
end
Cost(isnan(Cost)) = cost_thr;
Cost(:,:,layers+1) = single(~isnan(D));
CostF = recursive_bilateral(Cost, single(I), sigma_color, alpha_spacial);
clear Cost;
CostF = CostF(:,:,1:layers)./repmat(CostF(:,:,layers+1), [1 1 layers]);
D_hat = wta_simple(CostF, mindisp, 1);
clear CostF;

