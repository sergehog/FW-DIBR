clear
close all
clc

dataset = 'Art'; 
%dataset = 'Aloe'; 
%dataset = 'Reindeer'; 
%dataset = 'Moebius'; 
%dataset = 'Dolls'; 
%dataset = 'Teddy'; 
%dataset = 'Cones'; 


I1 = single(imread(['datasets\',dataset,'\view1.png']));
D1 = single(imread(['datasets\',dataset,'\disp1_noholes.png']));

[h, w] = size(D1);
mindisp = floor(min(D1(:)));
maxdisp = ceil(max(D1(:)));

if strcmp(dataset,'Teddy') ~= 0 || strcmp(dataset,'Cones') ~= 0        
    baseline = 163;
else
    baseline = 160; % mm        
end

f = 3740;

% Camera intrinsics matrix:
K = [f, 0, w/2;...
    0, f, h/2;...
    0, 0,   1];

%[Dx, Dy] = gradient(D1);
    
Z1 = f*baseline./(D1); %Z1 = f*baseline./(D1+200);
C1 = single(K*[eye(3), [0 0 0]']);

addpath(genpath('../matting'));

for i=1:6
    Is{i} = imread(['datasets\',dataset,'\view',num2str(i),'.png']);
end

%if ~exist(['saves\renders\',dataset,'_',num2str(2),'_true.png'], 'file')
%    for v=2:6    
%        offR = 50*(v-1);
%        Id = Is{v};
%        imwrite(uint8(Id(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_true.png'], 'png');
%    end
%end

%%
%figure; imshow(uint8(Is{5})); title('Desired Image');

%f = 3740;

baseline = 163;
Z1 = f*baseline./(D1);
K = [f, 0, (w-1)/2;...
    0, f, (h-1)/2;...
    0, 0,   1];

C1 = single(K*[eye(3), [0 0 0]']);
C5 = single(K*[eye(3), [-160 0 4]']);

[Iv, Zv, AlphaV, Df, Db, Ifg, Ibg1, Ibg, alpha, Dbg, Ibgv, Ifgv] = rendering_proposed2(C5, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);

%[Iv, ~, ~] = free_rendering(C5, I1, Z1, C1);
psnrA = 20*log10(255/sqrt(nanmean((Iv(:)-single(Is{5}(:))).^2)));
figure; imshow(uint8(Iv)); title(['Rendered Image, PSNR=', num2str(psnrA)]);

%%
ZNear = min(Z1(:));
ZFar = max(Z1(:));
D1 = floor(255*(ZNear./Z1) .* ((ZFar-Z1)/(ZFar-ZNear)) + 0.5);
%maxdisp = 255;
%mindisp = 0;
[Df, Db, Ifg, Ibg, alpha] = layered_decomposition(I1, D1, sigma_depth, errosion_depth, errosion_fg, errosion_bg);
figure; imshow(Df, [])
figure; imshow(Db, [])
figure; imshow(uint8(Ifg))
figure; imshow(alpha)
%% Proposed method
clear PSNRS SSIMS PSNRHVS
sigma_depth = 0.93;
hf_alpha = 0.67;
errosion_depth = 1;
errosion_fg = 1;
errosion_bg = 3;


a = figure;
for v=2:6    
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);
    tic
    %[I1out] = rendering_layered_inpainted(Cv, I1, D1, C1, f_baseline, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda);
    [I1out] = rendering_proposed(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);     
    toc
    %[I1out, Zv] = rendering_conventional(Cv, I1, Z1, C1, 0.8, 2);        
    psnr1 = psnr(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    ssim1 = ssim(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    [psnrhv1,~] = psnrhvsm_mex(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);

    PSNRS(v-1, 1) = psnr1;    
    SSIMS(v-1, 1) = ssim1;    
    PSNRHVS(v-1, 1) = psnrhv1;    
    figure(a);     
    imshow(uint8(I1out(:,1:end-offR,:))); title(['Proposed method ',num2str(v),', PSNR=', num2str(psnr1)]); drawnow;
    %imwrite(uint8(I1out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_proposed.png'], 'png');
    %figure; imshow(Id(:,1:end-offR,:))
end
%%

%% Conventional rendering + recursive compensation
clear PSNRS SSIMS PSNRHVS
a = figure;
for v=2:6    
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);
    [I1out, Zv] = rendering_conventional(Cv, I1, Z1, C1, 0.8, 2);        
    psnr1 = psnr(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    ssim1 = ssim(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    [psnrhv1,~] = psnrhvsm_mex(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);

    PSNRS(v-1, 1) = psnr1;    
    SSIMS(v-1, 1) = ssim1;    
    PSNRHVS(v-1, 1) = psnrhv1;    
    figure(a);     
    imshow(uint8(I1out(:,1:end-offR,:))); title(['Conventional ',num2str(v-1),', PSNR=', num2str(psnr1)]); drawnow;
    imwrite(uint8(I1out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_conventional.png'], 'png');
    %figure; imshow(Id(:,1:end-offR,:))
end
%% Layered Rendering (compensation as previously)
avg_alpha = 0.94;
hf_alpha = 0.7;
a = figure;

for v=2%2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);

    I2out = rendering_layered(Cv, I1, Z1, C1, avg_alpha, hf_alpha);
    psnr2 = psnr(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    ssim2 = ssim(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    [psnrhv2,~] = psnrhvsm_mex(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    PSNRS(v-1, 2) = psnr2;    
    SSIMS(v-1, 2) = ssim2;    
    PSNRHVS(v-1, 2) = psnrhv2;    
    figure(a); imshow(uint8(I2out(:,1:end-offR,:))); title(['Layered, ',num2str(v-1),'; PSNR=', num2str(psnr2)]); drawnow;
    imwrite(uint8(I2out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_layered.png'], 'png');
end
%% Layered Rendering + Matting
avg_alpha = 0.96;
hf_alpha = 0.7;
matting_radius = 1;
matting_sigma = 10;
matting_epsilon = 1e-4;
matting_lambda = 0.1;
errosion = 2;
a = figure;
for v=2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);
    [I3out] = rendering_layered_matting(Cv, I1, Z1, D1, C1, errosion, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda);
    %I2out = rendering_layered(Cv, I1, Z1, C1, avg_alpha, hf_alpha, 2);
    psnr3 = psnr(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    ssim3 = ssim(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    [psnrhv3,~] = psnrhvsm_mex(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    PSNRS(v-1, 3) = psnr3;    
    SSIMS(v-1, 3) = ssim3;    
    PSNRHVS(v-1, 3) = psnrhv3;    
    figure(a); imshow(uint8(I3out(:,1:end-offR,:))); title(['Layered, ',num2str(v-1),'; PSNR=', num2str(psnr3)]); drawnow;
    imwrite(uint8(I3out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_matting.png'], 'png');
end
%% Layered Rendering + Matting + Inpainting

%% Compensating Depth and Color images
r = 2; % errosion radius
%avg_alpha = 0.96;
%avg_alpha2 = 0.98;
t_thr = 1;
D_hat = recursive_gaussian(D1, avg_alpha)./recursive_gaussian(ones([h w], 'single'), avg_alpha);
FgA = D1 > (D_hat);
FgB = D1 > (D_hat + t_thr);
BgB = imerode(~FgB, strel('disk', r));

%figure; imshow(FgA); title('Foreground A');
%figure; imshow(BgB); title('Background B');

[Dx, Dy] = gradient(D1);
W = ~(D1 > D_hat);
W = single(imerode(W, strel('disk', r)));
%figure; imshow(W); title('Weight');

T2 = FgA + 1 - (~FgA);

figure; imshow(imerode(T2==2, strel('disk', r)) - imerode(T2==0, strel('disk', r)), [-1 1]); title('Trimap 2');

%r = 3;
filename_errosion = ['saves/',dataset,'_errosion', num2str(r),'.mat'];
if exist(filename_errosion, 'file')
    eval(['load ', filename_errosion]);
else
    disp(['errode# ', num2str(r)]);    
    Ti = uint8(imerode(T2==2, strel('disk', r)) + 1 - imerode(T2==0, strel('disk', r)));
    clear mex
    tic
    alpha = learning_matting(I1, Ti, 1, 1);
    toc
    
    [Fi, Bi] = solveFB(double(I1)/255, double(alpha));
    Fi = single(Fi*255);
    Bi = single(Bi*255);    
    
    % simple compensation
    F2i = Fi;
    F2i(repmat(alpha < 0.5, [1 1 3])) = nan;    
    F2i = recursive_compensation_simple(F2i, hf_alpha);
    %figure; imshow(uint8(Fi)); title('Given FG with holes');
    
    %B2i = Bi;
    %B2i(repmat(alpha > 0, [1 1 3])) = nan;
    B2i = I1;
    B2i(repmat(~BgB, [1 1 3])) = nan;
    B2i = recursive_compensation_simple(B2i, hf_alpha);
    figure; imshow(uint8(B2i)); title('Recursive-Compensated BG');

    % compensation with inpainting
    B3i = I1;
    %B3i = Bi;
    B3i(:,:,4) = Dx;
    B3i(:,:,5) = Dy;
    %B3i(repmat(alpha > 0, [1 1 5])) = nan;
    B3i(repmat(~BgB, [1 1 5])) = nan;
    
    W(W < 1e-4) = 1e-4;
    tic
    [B3i,~] = criminisi_inpainting(B3i, W, 0, 9, 'et', 'ag');
    toc    
    B3i = B3i(:,:,1:3);
    if sum(isnan(B3i(:))) > 0
        [B3i,~] = criminisi_inpainting(B3i, ones([h w], 'single'), 500, 5, 'et', 'ag');
    end
    figure; imshow(uint8(B3i(:,:,1:3))); title('Inpainted BG');       
    %figure; imshow(isnan(B3i(:,:,1)));
    % compensation of depth maps
    Df = D1;
    Df(T2 ~= 2) = 0;        
    Df = recursive_gaussian(Df, hf_alpha)./recursive_gaussian(single(T2 == 2), hf_alpha);
    %Df = recursive_compensation_simple(Df, hf_alpha);
    Zf = f*baseline./(Df);
    clear Df;
    
    Db = D1;
    Db(~BgB) = 0;        
    Db = recursive_gaussian(Db, hf_alpha)./recursive_gaussian(single(BgB), hf_alpha);  
    %Db = recursive_compensation_simple(Db, hf_alpha);        
    Db(T2 == 0) = D1(T2 == 0);    
    Zb = f*baseline./(Db);
    clear Db;
    
    eval(['save ', filename_errosion, ' Ti alpha Fi F2i Bi B2i B3i Zf Zb']);
    
    figure; imshow(uint8(F2i.*repmat(alpha, [1 1 3]) + repmat(reshape([255 150 50], [1 1 3]), [h w]).*repmat(1-alpha, [1 1 3]))); title(['Foreground with Alpha r=',num2str(r)]);
    drawnow;
    
end    
%%
a = figure;
for v=2:6
    disp(['view# ', num2str(v)]);
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);
        
    [Fv, ~, ~] = free_rendering(Cv, F2i, Zf, C1);
    [Av, ~, ~] = free_rendering(Cv, alpha, Zf, C1);    
    [B2v, Zbv, ~] = free_rendering(Cv, B2i, Zb, C1); % blurry
    [B3v, ~, ~] = free_rendering(Cv, B3i, Zb, C1); % inpainted    
    
    [B2v] = recursive_compensation(B2v, Zbv, hf_alpha);
    [B3v] = recursive_compensation(B3v, Zbv, hf_alpha);
    %figure; imshow(uint8(B3v))
    
    Fv(isnan(Fv)) = 0;
    Av(isnan(Av)) = 0;
    
    C2v = Fv.*repmat(Av, [1 1 3]) + B2v.*repmat((1-Av), [1 1 3]);        
    psnr3 = psnr(uint8(round(C2v(:,1:end-offR,:))), Id(:,1:end-offR,:));
    ssim3 = ssim(uint8(round(C2v(:,1:end-offR,:))), Id(:,1:end-offR,:));
    [psnrhv3,~] = psnrhvsm_mex(uint8(round(C2v(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);

    C3v = Fv.*repmat(Av, [1 1 3]) + B3v.*repmat((1-Av), [1 1 3]);        
    psnr4 = psnr(uint8(round(C3v(:,1:end-offR,:))), Id(:,1:end-offR,:));
    ssim4 = ssim(uint8(round(C3v(:,1:end-offR,:))), Id(:,1:end-offR,:));
    [psnrhv4,~] = psnrhvsm_mex(uint8(round(C3v(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    
    PSNRS(v-1, 3) = psnr3;    
    PSNRS(v-1, 4) = psnr4; 
    SSIMS(v-1, 3) = ssim3;    
    SSIMS(v-1, 4) = ssim4; 
    PSNRHVS(v-1, 3) = psnrhv3;    
    PSNRHVS(v-1, 4) = psnrhv4;    
    
    figure(a); imshow(uint8(C3v(:,1:end-offR,:))); title(['Inpainting + Matting, ',num2str(v-1),'; PSNR=', num2str(psnr4)]); drawnow;
    
    imwrite(uint8(C2v(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_matting.png'], 'png');
    imwrite(uint8(C3v(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_matting_inpainting.png'], 'png');
end    
    
%%
%PSNRS = PSNRS(:,1:4)
names = {'Conventional rendering', 'Layered Rendering', 'Layered + Matting', 'Inpainting + Matting'};
fig1 = figure; 
plot(PSNRS, 'LineWidth', 2); 
legend(names);
title(['Rendering performance for "',dataset,'" dataset']);
ylabel('PSNR, dB');
xlabel('Rendering offset');
%ylim([25 32]);
%ylim([22 33]);
saveas(fig1, ['figures/',dataset,'_psnr.pdf'], 'pdf');
saveas(fig1, ['figures/',dataset,'_psnr.fig'], 'fig');
%
fig2 = figure; 
plot(SSIMS, 'LineWidth', 2); 
legend(names);
title(['Rendering performance for "',dataset,'" dataset']);
ylabel('SSIM');
xlabel('Rendering offset');
saveas(fig2, ['figures/',dataset,'_ssim.pdf'], 'pdf');
saveas(fig2, ['figures/',dataset,'_ssim.fig'], 'fig');

%
fig3 = figure; 
plot(PSNRHVS, 'LineWidth', 2); 
legend(names);
title(['Rendering performance for "',dataset,'" dataset']);
ylabel('PSNR-HVS, dB');
xlabel('Rendering offset');
saveas(fig3, ['figures/',dataset,'_hvs.pdf'], 'pdf');
saveas(fig3, ['figures/',dataset,'_hvs.fig'], 'fig');

%ylim([25 31]);
%ylim([22 33]);
