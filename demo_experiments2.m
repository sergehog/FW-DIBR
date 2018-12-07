close all
clear 
clc

%dataset = 'Cones'; mindisp = 60; 
dataset = 'Art'; mindisp = 75; 
%dataset = 'Aloe'; mindisp = 40; 
%dataset = 'Reindeer'; mindisp = 60; 
%dataset = 'Moebius'; 


I1 = single(imread(['datasets\',dataset,'\view1.png']));
D1 = single(imread(['datasets\',dataset,'\disp1_noholes.png']));

[h, w] = size(D1);
%HW = h * w;
maxdisp = ceil(max(D1(:)));
D1(D1 < mindisp) = nan; % true depth with possible holes
%D1f = D1; % depth with compensated holes
%D1f(isnan(D1)) = 0;
%W = single(maxdisp-D1f);
%W(isnan(D1)) = 0;
%D1f = hypothesis_filter(D1, I1, 0.1, 0.9);
%D1f = recursive_gaussian(D1f.*W, 0.8)./recursive_gaussian(W, 0.8);
%D1f(~isnan(D1)) = D1(~isnan(D1));

f = 3740;
% Camera intrinsics matrix:
K = [f, 0, (w)/2;...
    0, f, (h)/2;...
    0, 0,   1];
baseline = 160; % mm    
%[Dx, Dy] = gradient(D1);
    
Z1 = f*baseline./(D1); %Z1 = f*baseline./(D1+200);
%Z1f = f*baseline./(D1f);
C1 = single(K*[eye(3), [0 0 0]']);
Zmin = min(Z1(:));
Zmax = max(Z1(:));


for i=1:6
    Is{i} = imread(['datasets\',dataset,'\view',num2str(i),'.png']);
end

if ~exist(['saves\renders\',dataset,'_',num2str(2),'_true.png'], 'file')
    for v=2:6    
        offR = 50*(v-1);
        Id = Is{v};
         imwrite(uint8(Id(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_true.png'], 'png');
    end
end
addpath('D:\Work\2015\paperEI\MVQMuX');
%figure; imshow(uint8(I1)); title('Given Image');
figure; imshow(D1, [mindisp maxdisp]); colormap(pink); title('Given Depth');
%figure; imshow(D1f, [mindisp maxdisp]); colormap(pink); title('Hole-filled Depth');

v=6; % desired view index
offR = 50*(v-1);
Id = Is{v};
t = [-40*(v-1) 0 0]';
Cv = single(K*[eye(3), t]);    
figure; imshow(uint8(Id(:,1:end-offR,:))); title('Desired Image');
%%
[Iv, Zv, ~] = free_rendering(Cv, I1, Z1, C1);
if radius > 0
    Valid = imerode(~isnan(Zv),strel('disk',radius));
    Iv(repmat(~Valid, [1 1 3])) = nan;
    Zv(~Valid) = nan;
end
figure; imshow(uint8(Iv(:,1:end-offR,:))); title('Rendered with no filling');
imwrite(uint8(Iv(:,1:end-offR,:)), 'saves/renderd_no_filling.png', 'png');
%%
sigma_depth = 0.92;
hf_alpha = 0.67;
errosion_depth = 1;
errosion_fg = 1;
errosion_bg = 3;

[Df, Db, Ifg, Ibg, alpha] = layered_decomposition(I1, D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg);


figure; imshow(uint8(Ifg)); title('Ifg');
figure; imshow(uint8(Ibg)); title('Ibg');
figure; imshow(alpha, [0 1]); title('Alpha');

%figure; imshow(Db, [0 255]); colormap(pink); title('Db');
%figure; imshow(Df, [0 255]); colormap(pink); title('Df');

%%

figure; imshow(alpha, []); title('Estimated Alpha');
    
Dbg = recursive_compensation_weighted(Db, (255-Db), hf_alpha);
Dbg(Dbg > D1) = D1(Dbg > D1)-1;
figure; imshow(Dbg, [0 255]); colormap(pink); title('Compensated Bg Depth');

Pr = (max(Dbg(:))-Dbg)/max(Dbg(:));
Pr(Pr<1e-5) = 1e-5;
Pr(isnan(Pr)) = 1e-5;
Pr(~isnan(Df)) = 1e-5;

%figure; imshow(Pr, [])
%%
tic
[Ib2] = inpaint_rgbd(Ibg, Dbg, Pr, 300, 0);
toc
figure; imshow(uint8(Ib2))
%figure; imshow(Ibg(:,:,1), [])

%% PROPOSED
clc
sigma_depth = 0.93;
hf_alpha = 0.67;
errosion_depth = 1;
errosion_fg = 1;
errosion_bg = 3;

tic
[Iv, Zv, Av] = rendering_proposed(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);
toc
psnr1 = psnr(uint8(round(Iv(:,1:end-offR,:))), uint8(Id(:,1:end-offR,:)));
figure; imshow(uint8(Iv(:,1:end-offR,:))); title(['Rendered with proposed approach, PSNR=', num2str(psnr1)]);
%figure; imshow(Av(:,1:end-offR,:), [])
%% Layered Approach (Simple Hole-Filling)
%close all

sigma_depth = 0.93;
hf_alpha = 0.67;
errosion_depth = 1;
errosion_fg = 1;
errosion_bg = 3;

[Ib, Zb, Alpha] = rendering_layered(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);
psnrB = psnr(uint8(round(Ib(:,1:end-offR,:))), Id(:,1:end-offR,:));
figure; imshow(uint8(Ib(:,1:end-offR,:))); title(['Layered Rendering, PSNR=', num2str(psnrB)]); drawnow;
figure; imshow(Zb, []); colormap(pink);
figure; imshow(Alpha(:,1:end-offR,:), [0 1]); title('Alpha');
%figure; imshow(isnan(Ib(:,:,1)))


    
%% Rendering with Hierarchical Hole-filling
[Ia, Za] = rendering_hierarchical(Cv, I1, Z1, C1, 3);
psnrA = psnr(uint8(round(Ia(:,1:end-offR,:))), Id(:,1:end-offR,:));
figure; imshow(uint8(Ia(:,1:end-offR,:))); title(['Hierarchical Hole-Filling, PSNR=', num2str(psnrA)]); drawnow;

%% JTDI
search = 100;
patch = 4;
errorsion = 2;
tic
[Iv, Dv, Conf, Data] = rendering_jtdi(Cv, I1, Z1, C1, search, patch, errorsion);
toc
psnr1 = psnr(uint8(round(Iv(:,1:end-offR,:))), uint8(Id(:,1:end-offR,:)));
clear mex
figure; imshow(uint8(Iv(:,1:end-offR,:))); title(['Rendered with JTDI, PSNR=', num2str(psnr1)]);

%% Rendering DE WITH


patch = 7;
search = 200;
errorsion = 2;
fg_radius = 100;
tic

Iv = rendering_dewith(Cv, I1, Z1, C1, fg_radius, patch, search, errorsion);
toc
clear mex
psnr1 = psnr(uint8(round(Iv(:,1:end-offR,:))), uint8(Id(:,1:end-offR,:)));
figure; imshow(uint8(Iv(:,1:end-offR,:))); title(['Rendered with De With approach, PSNR=', num2str(psnr1)]);












%%
errosion_depth = 0;
avg_alpha = 0.96;
thr0 = 0;
thr1 = 1;

[h, w] = size(D1);
HW = h * w;
D = D1;
D(isnan(D)) = 0;
D_hat = recursive_gaussian(D, avg_alpha)./recursive_gaussian(single(~isnan(D1)), avg_alpha);
Mask0 = D1 > (D_hat + thr0) & ~isnan(D1);
Mask1 = D1 > (D_hat + thr1) & ~isnan(D1);
MDiff = Mask0 & ~Mask1;

FgDepthMask = imerode(Mask0 | MDiff, strel('disk', errosion_depth)) | MDiff;
BgDepthMask = imerode(~Mask1 | MDiff, strel('disk', errosion_depth)) | MDiff;

%FgA = D1 > (D_hat + t_thr) & ~isnan(D1);

% some overlap between Bg and Fg is required
%Both = imopen( (D1 > D_hat)  & ~(D1 > D_hat + 3), strel('disk', t_thr));
%Both = (D1 > D_hat)  & ~(D1 > D_hat + 3) & ~isnan(D1);

%FgDepthMask = FgA | Both;
%BgDepthMask = imerode(~(D1 > D_hat), strel('disk', 1)) | Both;
%BgDepthMask = imerode(~(D1 > D_hat) & ~isnan(D1)) | Both;
%BgDepthMask = imerode((~(D1 > D_hat) & ~isnan(D1)) | Both, strel('disk', 1)) | Both;
%figure; imshow(FgDepthMask); title('FgDepthMask');
%figure; imshow(BgDepthMask); title('BgDepthMask');
D1f = D1;
D1f(~FgDepthMask) = nan;
Db = D1;
Db(~BgDepthMask) = nan;
figure; imshow(D1f, [mindisp maxdisp]); colormap(pink); title('Fg Depth');
figure; imshow(Db, [mindisp maxdisp]); colormap(pink); title('Bg Depth');
%%
matting_radius = 2;
matting_sigma = 10;
matting_epsilon = 1e-4;
matting_lambda = 0.01;
errosion_fg = 1;
errosion_bg = 2;
alpha_thr = 0.01;
FgColorMask = imerode(Mask0 & ~isnan(D1), strel('disk', errosion_fg)) | MDiff;
BgColorMask = (imerode((~Mask1 & ~isnan(D1)), strel('disk', errosion_bg)) | MDiff) & (~FgColorMask);
Ti = FgColorMask + 1 - BgColorMask;
figure; imshow(Ti, [0 2]); title('Color Trimap');
%Bg = I1;
%Bg(repmat(~BgColorMask & ~Both, [1 1 3])) = nan;
%figure; imshow(uint8(Bg))
%%
    AL = get_learning_laplacian(uint8(I1), uint8(Ti), matting_radius, matting_sigma, matting_epsilon);
    S = double(Ti(:)~=1);
    S = spdiags(S(:), 0, HW, HW);
    alpha = (matting_lambda*AL + S'*S) \ (S'*(double(Ti(:))));
    alpha = single(reshape(alpha, [h w]))/2;
    %alpha = single(reshape(alpha, [h w]) + 1)/2;
    alpha(alpha < 0) = 0;
    alpha(alpha > 1) = 1;
    figure; imshow(alpha, [0 1]); title(['Reconstructed Alpha, errosion-fg=', num2str(errosion_fg), '; errosion-bg=', num2str(errosion_bg)]);
%%
    [Fi, Bi] = solveFB(double(I1)/255, double(alpha));
    Fi = single(Fi*255);
    Bi = single(Bi*255);    
    Fi(repmat(MDiff, [1 1 3])) = I1(repmat(MDiff, [1 1 3]));
    Bi(repmat(MDiff, [1 1 3])) = I1(repmat(MDiff, [1 1 3]));
    
    %figure; imshow(uint8(Fi)); title('Foreground');
    figure; imshow(uint8(repmat(alpha, [1 1 3]) .* Fi + repmat(1-alpha, [1 1 3]) * 127)); title('Foreground');
    %figure; imshow(uint8(Bi)); title('Background');
    Bg = Bi;
    Bg(repmat(BgDepthMask, [1 1 3])) = I1(repmat(BgDepthMask, [1 1 3]));
    figure; imshow(uint8(Bg)); title('Background 2');
    Bg(repmat(~BgDepthMask, [1 1 3])) = nan;
    Bg(repmat(isnan(Db), [1 1 3])) = nan;
    figure; imshow(uint8(Bg)); title('Known Background');
%%
    Pr = (maxdisp-Db)/(maxdisp-mindisp);
    Pr(MDiff) = 1e-5;
    Pr(isnan(Pr)) = 1e-5;
    %[Dx, Dy] = gradient(D1);
    Bi = Bg;
    Bi(:,:,4) = D;
    %Bi(:,:,5) = Dy;
    Bi(repmat(~BgDepthMask, [1 1 4])) = nan;
    %B3i(repmat(alpha > 0, [1 1 5])) = nan;
    %B3i(repmat(~BgB, [1 1 5])) = nan;    
    %W(W < 1e-4) = 1e-4;
    
    tic
    [Bgi,~] = criminisi_inpainting(Bi, Pr, 0, 8, 'et', 'ag');
    toc    
    figure; imshow(uint8(Bgi(:,:,1:3))); title('Criminisi inpainting');
    figure; imshow(Bgi(:,:,4), [mindisp maxdisp]); colormap(pink); title('Inpainted Depth');
    %figure; imshow(Bgi(:,:,4), [-1 1]); colormap(jet); title('Inpainting Dx');
    %figure; imshow(Bgi(:,:,5), [-1 1]); colormap(jet); title('Inpainting Dy');
%%
tic
[alpha, If, Ib, D1f, Db] = prepare_matting_inpainting(I1, D1, errosion, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda);
toc
figure; imshow(alpha); title('Estimated Alpha');
figure; imshow(uint8(repmat(alpha,[1 1 3]).*If)); title('Foreground Image');
figure; imshow(uint8(Ib)); title('Background Image');
figure; imshow(D1f, [mindisp maxdisp]); colormap(pink); title('Foreground Depth');
figure; imshow(Db, [mindisp maxdisp]); colormap(pink); title('Background Depth');
%Db(isnan(Ib(:,:,1))) = nan;
%%
    %alpha = 0.8;
    
    %[Iout, Ar, Ifr, Ibr, If, Ib, Zf, Zb, alpha] = rendering_layered_matting(Cv, I1, Z1, D1, C1, errosion, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda);
    %[Iout, Zv] = rendering_conventional(Cv, I1, Z1, C1, hf_alpha, errosion);    
    %[Iout] = rendering_layered(Cv, I1, Z1, C1, avg_alpha, hf_alpha, 2);    
    ypsnr1 = psnr(rgb2gray(uint8(round(Iout(:,1:end-offR,:)))), rgb2gray(uint8(Id(:,1:end-offR,:))));
    
    psnr1 = psnr(uint8(round(Iout(:,1:end-offR,:))), uint8(Id(:,1:end-offR,:)));
    %mp_psnr = MP_PSNR_reduc(uint8(Iout(:,1:end-offR,:)), uint8(I1(:,1:end-offR,:)));
    [mp_psnr, dswim] = MVQMuX(uint8(Iout(:,1:end-offR,:)), uint8(I1(:,1:end-offR,:)));
    %mp_psnr = MP_PSNR_reduc(uint8(I1(:,1:end-offR,:)), uint8(Iout(:,1:end-offR,:)));
    %psnr1 = 0;
    figure; imshow(uint8(Iout(:,1:end-offR,:))); title(['Rendered View; PSNR=', num2str(psnr1), '; MP-PSNR=', num2str(mp_psnr),'; 3DSWIM=', num2str(dswim)]);
    %figure; imshow(Zv, [Zmin Zmax]); colormap(pink); title('Rendered Depth');
    %figure; imshow(Ar(:,1:end-offR), [0 1])
    %figure; imshow(uint8(If))
    %figure; imshow(uint8(Ib))
    %figure; imshow(Zb, []); colormap(jet)
    %figure; imshow(Zf, []); colormap(jet)
    
