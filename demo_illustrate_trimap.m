close all
clear 
clc

I = imread('D:\Work\stereo\Art\imL.png');
D = single(imread('D:\Work\stereo\Art\disp.png'));

mindisp = 70;
D (D < mindisp) = nan;
maxdisp = max(D(:));
figure; imshow(I);
figure; imshow(D, [mindisp maxdisp]);% colormap(gca, pink);
%%
beta_R = 60;
beta_F = 0.2;

Davg = fast_average(D, beta_R);
%figure; imshow(Davg, [mindisp maxdisp]);% colormap(gca, pink);
%title('Averaged Depth');
F = D > (Davg + beta_F);
figure; imshow(F); title('Fg');

B = 1 - F;
figure; imshow(B); title('Bg');

%%
beta_K = 2;
Fs = imerode(F, strel('disk', beta_K));
Bs = imerode(B, strel('disk', beta_K));
Ts = Fs - Bs + 1;
%%
figure; imshow(Ts(300:1100,700:1300), [0 2]); title('Trimap');
imwrite(uint8(round(0.5*Ts*255)), 'saves/trimap_normal.png', 'png');
Ts2 = Ts(300:1100,700:1300);
imwrite(uint8(round(0.5*Ts2*255)), 'saves/trimap2_normal.png', 'png');
%%
beta_G = 10;
beta_R = 1;
[Dx, Dy] = gradient(D);
Uini = Dx.^2 + Dy.^2 > beta_G;
U = imdilate(Uini, strel('disk', beta_R));
figure; imshow(U, [0 1]); title('Unknown areas');
imwrite(U, 'saves/trimap_unknown.png', 'png');
U2 = U(300:1100,700:1300);
imwrite(U2, 'saves/trimap2_unknown.png', 'png');

%%

Ff = F - (F-Fs).*U;
Bf = B - (B-Bs).*U;
Tf = Ff - Bf + 1;
figure; imshow(Tf, [0 2]); title('New Trimap');

%Fs = imerode(F, strel('disk', beta_K));
%Bs = imerode(B, strel('disk', beta_K));
%Ts = Fs - Bs + 1;
%figure; imshow(Ts, [0 2]); title('Trimap');
%figure; imshow(Tf(300:1100,700:1300), [0 2]); title('New Trimap');
%imwrite(uint8(round(0.5*Tf*255)), 'saves/trimap_new.png', 'png');
%Tf2 = Tf(300:1100,700:1300);
%imwrite(uint8(round(0.5*Tf2*255)), 'saves/trimap2_new.png', 'png');

%% 
%L = get_laplacian(single(I), 1, 20, 1e-4);
L = get_learning_laplacian(I, uint8(Ts), 1, 20, 1e-4);
[h, w, ~] = size(I);
HW = h * w;
%%

y = double(Tf(:)-1);
S = spdiags(double(Tf(:)~=1), 0, HW, HW);
lambda = 0.001;
tic
x = (S + lambda*(L)) \ (y);
toc
figure; imshow(reshape(x, [h w]), [-1 1]); title('Alpha map');
%%
alpha = (reshape(x, [h w]) + 1)/2;
alpha(alpha < 0) = 0;
alpha(alpha > 1) = 1;
figure; imshow(alpha(300:1100,700:1300), [0 1]); title('Alpha map');
imwrite(uint8(round(alpha*255)), 'saves/alpha_new.jpg', 'jpg');
alpha2 = alpha(300:1100,700:1300);
imwrite(uint8(round(alpha2*255)), 'saves/alpha2_new.jpg', 'jpg');
%%

[F, B] = solveFB(im2double(I), alpha);
%F1 = I.*repmat(alpha, [1 1 3]);
F2 = F.*repmat(alpha, [1 1 3]);

figure; imshow(F2)

imwrite(uint8(round(F2*255)), 'saves/fg_new.jpg', 'jpg');
alpha2 = alpha(300:1100,700:1300);
imwrite(uint8(round(alpha2*255)), 'saves/alpha2_new.jpg', 'jpg');
%%
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

[Iv, Zv, AlphaV, Df, Db, Ifg, Ibg1, Ibg, alpha, Dbg, Ibgv, Ifgv, AlphaV] = rendering_proposed(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);
