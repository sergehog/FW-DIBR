function [Iv, Zv, AlphaV] = rendering_proposed(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg, inpaint_radius, inpaint_search)
if nargin < 5
    sigma_depth = 0.9;    
end
if nargin < 6
    hf_alpha = 0.8;    
end
if nargin < 7
    errosion_depth = 0;
end
if nargin < 8
    errosion_fg = 1;    
end
if nargin < 9
    errosion_bg = 3;    
end
if nargin < 10
    inpaint_radius = 7;
end
if nargin < 11
    inpaint_search = 300;
end

ZNear = min(Z1(:));
ZFar = max(Z1(:));
D1 = floor(255*(ZNear./Z1) .* ((ZFar-Z1)/(ZFar-ZNear)) + 0.5);
%maxdisp = 255;
%mindisp = 0;
[Df, Db, Ifg, Ibg, alpha] = layered_decomposition(I1, D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg);
%[Df, Db, Ifg, Ibg, alpha] = layered_decomposition(D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg);
%figure; imshow(Df, [0 255]); colormap(pink); title('Fg Depth'); drawnow;
%figure; imshow(Db, [0 255]); colormap(pink); title('Bg Depth');drawnow;
%figure; imshow(uint8(Ifg)); title('Fg Color'); drawnow;
%figure; imshow(uint8(Ibg)); title('Bg Color'); drawnow;
%figure; imshow(alpha); title('alpha'); drawnow;

Ifga = recursive_compensation_simple(Ifg, hf_alpha);
Dfa = recursive_compensation_simple(Df, hf_alpha);
%figure; imshow(uint8(Ifga)); title('Compensated Fg Color'); drawnow;
%figure; imshow(Dfa, [0 255]); colormap(pink); title('Compensated Fg Depth'); drawnow;
Zfg = 1./((Dfa/255)*(1/ZNear - 1/ZFar) + 1/ZFar);   
clear Df Dfa;
%return
%Zfg = 1./((Df/255)*(1/ZNear - 1/ZFar) + 1/ZFar);   
%Zfg = f_baseline./Df;
%clear Df;
alpha(isnan(alpha)) = 0;
alpha(alpha < 0) = 0;
alpha(alpha > 1) = 1;
%figure; imshow(alpha, [0 1]); title('alpha');drawnow;

%%
%alpha_thr = 0.01;
%FgColorMask = imerode(Mask0 & ~isnan(D1), strel('disk', errosion_fg)) | MDiff;
%BgColorMask = (imerode((~Mask1 & ~isnan(D1)), strel('disk', errosion_bg)) | MDiff) & (~FgColorMask);

%Ti = FgColorMask + 1 - BgColorMask;
%figure; imshow(Ti, []); title('Color Trimap'); drawnow;
%Bg = I1;
%Bg(repmat(~BgColorMask & ~Both, [1 1 3])) = nan;
%figure; imshow(uint8(Bg))
%%
%    AL = get_learning_laplacian(uint8(I1), uint8(Ti), matting_radius, matting_sigma, matting_epsilon);
%    S = double(Ti(:)~=1);
%    S = spdiags(S(:), 0, HW, HW);
%    alpha = (matting_lambda*AL + S'*S) \ (S'*(double(Ti(:))));
%    alpha = single(reshape(alpha, [h w]))/2;
%    %alpha = single(reshape(alpha, [h w]) + 1)/2;
%    %figure; imshow(alpha, [0 1]); title(['Reconstructed Alpha, errosion-fg=', num2str(errosion_fg), '; errosion-bg=', num2str(errosion_bg)]); drawnow;

%%
%    [Ifg, Bi] = solveFB(double(I1)/255, double(alpha));
%    Ifg = single(Ifg*255);
%    Bi = single(Bi*255);    
%    Ifg(repmat(alpha<=0.1, [1 1 3])) = nan;
%    Ifg(repmat(MDiff, [1 1 3])) = I1(repmat(MDiff, [1 1 3]));    
%    Ifg_f = Ifg;
%    Ifg_f(isnan(Ifg)) = 0;
%    Ifg_f = recursive_gaussian(Ifg_f, 0.8)./recursive_gaussian(single(~isnan(Ifg)), 0.8);
%    Ifg(isnan(Ifg)) = Ifg_f(isnan(Ifg));
%    clear Ifg_f;
%    Bi(repmat(MDiff, [1 1 3])) = I1(repmat(MDiff, [1 1 3]));
    
    %figure; imshow(uint8(Ifg)); title('Compensated Foreground');drawnow;
    %figure; imshow(uint8(repmat(alpha, [1 1 3]) .* Ifg + repmat(1-alpha, [1 1 3]) * 127)); title('Foreground'); drawnow;
    %figure; imshow(uint8(Bi)); title('Background');
%    Bg = I1;
%    Bg(repmat(~BgColorMask & ~MDiff, [1 1 3])) = nan;
    %Bg = Bi;
    %Bg(repmat(BgDepthMask, [1 1 3])) = I1(repmat(BgDepthMask, [1 1 3]));
    %figure; imshow(uint8(Bg)); title('Background 2');
    %Bg(repmat(~BgDepthMask, [1 1 3])) = nan;
%    Bg(repmat(isnan(Db), [1 1 3])) = nan;
    %figure; imshow(uint8(Bg)); title('Known Background'); drawnow;
%%
    
        
    Dbg = recursive_compensation_weighted(Db, (255-Db), hf_alpha);
    Dbg(Dbg > D1) = D1(Dbg > D1)-1;

    Pr = (max(Dbg(:))-Dbg)/(max(Dbg(:)-min(Dbg(:))));
    Pr(Pr<1e-5) = 1e-5;
    %Pr(MDiff) = 1e-5;
    Pr(isnan(Pr)) = 1e-5;
    
    %[Dx, Dy] = gradient(D1);
    %Bi = Ibg;
    %Bi(:,:,4) = Dba;
    %Bi(:,:,5) = Dy;
    %Bi(repmat(~BgDepthMask, [1 1 4])) = nan;
    %B3i(repmat(alpha > 0, [1 1 5])) = nan;
    %B3i(repmat(~BgB, [1 1 5])) = nan;    
    %W(W < 1e-4) = 1e-4;
    

    tic
    %[Bgi,~] = criminisi_inpainting(Bi, Pr, inpaint_search, inpaint_radius, 'et', 'ag');
    %[Bgi,~] = criminisi_inpainting(Bi, Pr, inpaint_search, inpaint_radius);
    %[Ibg, Dbg, ~, ~] = inpaint_rgbd_jtdi(Ibg, Dba, inpaint_search, inpaint_radius);
    [Ibg] = inpaint_rgbd(Ibg, Dbg, Pr, inpaint_search, inpaint_radius);    
    toc    
    %Ibg = Bgi(:,:,1:3);
    %Dbg = Bgi(:,:,4);
    %Zbg = f_baseline./Dbg;
    Dbg(Dbg > D1) = D1(Dbg > D1) - 1;
    Zbg = 1./((Dbg/255)*(1/ZNear - 1/ZFar) + 1/ZFar);   
    %figure; imshow(uint8(Ibg)); title('Inpainted Image'); drawnow;
    %figure; imshow(Dbg, [0 255]); colormap(pink); title('Inpainted Depth'); drawnow;
    clear Dbg Bgi Bi;
%%
[Ibgv, Zbv, ~] = free_rendering(Cv, Ibg, Zbg, C1);
[Ifgv, Zfv, ~] = free_rendering(Cv, Ifga, Zfg, C1);
[AlphaV, ~, ~] = free_rendering(Cv, single(alpha), Zfg, C1);

Ibgv = recursive_compensation_weighted(Ibgv, repmat(Zbv, [1 1 3]), hf_alpha);
Ifgv = recursive_compensation_simple(Ifgv, hf_alpha);
AlphaV = recursive_compensation_simple(AlphaV, hf_alpha);

Iv = Ibgv.*(1-repmat(AlphaV,[1 1 3])) + Ifgv.*repmat(AlphaV,[1 1 3]);

Zv = Zbv;
Zv(AlphaV > 0.3) = Zfv(AlphaV > 0.3);
 