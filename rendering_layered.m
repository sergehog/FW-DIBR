function [Iv, Zv, AlphaV] = rendering_layered(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg)
if nargin < 5
    sigma_depth = 0.88;
end
if nargin < 6
    hf_alpha = 0.88;
end
if nargin < 7
    errosion_depth = 1;
end
if nargin < 8
    errosion_fg = 1;
end
if nargin < 9
    errosion_bg = 2;
end

ZNear = min(Z1(:));
ZFar = max(Z1(:));
D1 = floor(255*(ZNear./Z1) .* ((ZFar-Z1)/(ZFar-ZNear)) + 0.5);

[Df, Db, Ifg, Ibg, Alpha] = layered_decomposition(I1, D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg);
%figure; imshow(Dba, [0 255]); colormap(pink); title('Dba');

Dba = recursive_compensation_weighted(Db, (255-Db), hf_alpha); 
Dba(Dba > D1) = D1(Dba > D1);
Ibga = recursive_compensation_weighted(Ibg, repmat(255-Db,[1 1 3]), hf_alpha);
Zbg = 1./((Dba/255)*(1/ZNear - 1/ZFar) + 1/ZFar);   
%figure; imshow(uint8(Ibga)); title('Ibga');
%figure; imshow(Dba, [0 255]); colormap(pink); title('Dba');
clear Db Dba;

Ifga = recursive_compensation_simple(Ifg, hf_alpha);
Dfa = recursive_compensation_simple(Df, hf_alpha);
Zfg = 1./((Dfa/255)*(1/ZNear - 1/ZFar) + 1/ZFar);   
%figure; imshow(uint8(Ifga)); title('Ifga');
%figure; imshow(Dfa, [0 255]); colormap(pink); title('Dfa');
clear Df Dfa;

%figure; imshow(Alpha, [0 1]); title('Alpha');


[Ibgv, Zbv, ~] = free_rendering(Cv, Ibga, Zbg, C1);
[Ifgv, Zfv, ~] = free_rendering(Cv, Ifga, Zfg, C1);
[AlphaV, ~, ~] = free_rendering(Cv, single(Alpha), Zfg, C1);

Ibgv = recursive_compensation_weighted(Ibgv, repmat(Zbv, [1 1 3]), hf_alpha);
Ifgv = recursive_compensation_simple(Ifgv, hf_alpha);
AlphaV = recursive_compensation_simple(AlphaV, hf_alpha);
%figure; imshow(uint8(Ifgv)); title('Rendered Fg');
%figure; imshow(Zfgv, []); title('Rendered Fg Z'); colormap(jet)
%Fg(isnan(Fg)) = 0;
%Fg(isnan(Ifgv(:,:,1))) = 0;
%FgMask = imclose(Fg>0, strel('disk', 1));
%FgMask = imopen(Fg>0, strel('disk', 2));
%Fg(~FgMask) = 0;
%figure; imshow(Fg, [0 1]); title('Rendered Alpha');


%Ifgv(isnan(Ifgv)) = 0;

%k = [0.5 0.8 0.5; 0.8 1 0.8; 0.5 0.8 0.5;];
%Fga = imfilter(Fg, k/sum(k(:)));
%Fga(Fg==0) = 0;
%AlphaV(AlphaV) = 0;

%figure; imshow(uint8(Ibgv)); title('Rendered & Compensated Background'); 
%figure; imshow(uint8(Ifgv)); title('Rendered & Compensated Foreground'); 

%BgHoles = isnan(Zbv) | isnan(Ibgv(:,:,1));
%Ibgv(repmat(BgHoles, [1 1 3])) = nan;
%Zbgv(BgHoles) = nan;
%[Ibgva, ~] = recursive_compensation(Ibgv, Zbgv, hf_alpha);
Iv = Ibgv.*(1-repmat(AlphaV,[1 1 3])) + Ifgv.*repmat(AlphaV,[1 1 3]);
Zv = Zbv;
Zv(AlphaV > 0.3) = Zfv(AlphaV > 0.3);

%figure; imshow(Fg, []); title('Rendered Foreground Mask');  colormap(pink)