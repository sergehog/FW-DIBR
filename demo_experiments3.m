clear
close all
clc

%dataset = 'Art'; 
%dataset = 'Reindeer'; 
%dataset = 'Aloe'; 
%dataset = 'Moebius'; 
dataset = 'Cones'; 
%dataset = 'Teddy'; 

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
%baseline = 160; % mm    
%[Dx, Dy] = gradient(D1);
    
Z1 = f*baseline./(D1); %Z1 = f*baseline./(D1+200);
C1 = single(K*[eye(3), [0 0 0]']);

addpath(genpath('../matting'));

for i=1:6
    Is{i} = imread(['datasets\',dataset,'\view',num2str(i),'.png']);
end

if ~exist(['saves\renders\',dataset,'_',num2str(2),'_true.png'], 'file')
    for v=2:6    
        offR = 50*(v-1);
        Id = Is{v};
        imwrite(uint8(Id(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_true.png'], 'png');
        save(['saves\mat\',dataset,'_',num2str(v),'_true.mat'], 'Id');
    end
end
%% Proposed method (layered with inpainting) 
% 32.6806
%clear PSNRS SSIMS PSNRHVS
sigma_depth = 0.93;
hf_alpha = 0.67;
errosion_depth = 1;
errosion_fg = 1;
errosion_bg = 3;

if 1==1
a = figure;
for v=2:6    
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);
    tic
    %[I1out] = rendering_layered_inpainted(Cv, I1, D1, C1, f_baseline, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda);
    [I1out] = rendering_proposed(Cv, I1, Z1, C1, sigma_depth, hf_alpha, errosion_depth, errosion_fg, errosion_bg);
    I1out = round(I1out);
    toc
    %[I1out, Zv] = rendering_conventional(Cv, I1, Z1, C1, 0.8, 2);        
    psnr1 = psnr(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    %ssim1 = ssim(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:));
    %[psnrhv1,~] = psnrhvsm_mex(uint8(round(I1out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);

    %PSNRS(v-1, 1) = psnr1;    
    %SSIMS(v-1, 1) = ssim1;    
    %PSNRHVS(v-1, 1) = psnrhv1;    
    figure(a);     
    imshow(uint8(I1out(:,1:end-offR,:))); title(['Proposed method ',num2str(v),', PSNR=', num2str(psnr1)]); drawnow;
    imwrite(uint8(I1out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_proposed.png'], 'png');
    save(['saves\mat\',dataset,'_',num2str(v),'_proposed.mat'], 'I1out');
    %figure; imshow(Id(:,1:end-offR,:))
end
end
%% Proposed method (layered with recursive compensation)
avg_alpha = 0.94;
hf_alpha = 0.7;
a = figure;

for v=2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);

    I2out = rendering_layered(Cv, I1, Z1, C1, avg_alpha, hf_alpha);
    I2out = round(I2out);
    psnr2 = psnr(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    %ssim2 = ssim(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    %[psnrhv2,~] = psnrhvsm_mex(uint8(round(I2out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    %PSNRS(v-1, 2) = psnr2;    
    %SSIMS(v-1, 2) = ssim2;    
    %PSNRHVS(v-1, 2) = psnrhv2;    
    figure(a); imshow(uint8(I2out(:,1:end-offR,:))); title(['Layered, ',num2str(v-1),'; PSNR=', num2str(psnr2)]); drawnow;
    imwrite(uint8(I2out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_layered.png'], 'png');
    save(['saves\mat\',dataset,'_',num2str(v),'_layered.mat'], 'I2out');
end
%% Conventional method (with HHF compensation)
a = figure;

for v=2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);

    [I3out, ~] = rendering_hierarchical(Cv, I1, Z1, C1, 3);
    I3out = round(I3out);
    psnr3 = psnr(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    %ssim3 = ssim(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    %[psnrhv3,~] = psnrhvsm_mex(uint8(round(I3out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    %PSNRS(v-1, 3) = psnr3;    
    %SSIMS(v-1, 3) = ssim3;    
    %PSNRHVS(v-1, 3) = psnrhv3;    
    figure(a); imshow(uint8(I3out(:,1:end-offR,:))); title(['HHF, ',num2str(v-1),'; PSNR=', num2str(psnr3)]); drawnow;
    imwrite(uint8(I3out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_hhf.png'], 'png');
    save(['saves\mat\',dataset,'_',num2str(v),'_hhf.mat'], 'I3out');
end
%% JTDI
search = 300;
patch = 4;
errorsion = 2;
a = figure;
for v=2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);

    [I4out, ~, ~, ~] = rendering_jtdi(Cv, I1, Z1, C1, search, patch, errorsion);    
    I4out = round(I4out);
    psnr4 = psnr(uint8(round(I4out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    %ssim4 = ssim(uint8(round(I4out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    %[psnrhv4,~] = psnrhvsm_mex(uint8(round(I4out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    %PSNRS(v-1, 4) = psnr4;    
    %SSIMS(v-1, 4) = ssim4;    
    %PSNRHVS(v-1, 4) = psnrhv4;    
    figure(a); imshow(uint8(I4out(:,1:end-offR,:))); title(['JTDI, ',num2str(v-1),'; PSNR=', num2str(psnr4)]); drawnow;
    imwrite(uint8(I4out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_jtdi.png'], 'png');   
    save(['saves\mat\',dataset,'_',num2str(v),'_jtdi.mat'], 'I4out');
end

%% DeWith Method
patch = 7;
search = 200;
errorsion = 2;
fg_radius = 100;
a = figure;
for v=2:6        
    offR = 50*(v-1);
    Id = Is{v};
    t = [-40*(v-1) 0 0]';
    Cv = single(K*[eye(3), t]);

    I5out = rendering_dewith(Cv, I1, Z1, C1, fg_radius, patch, search, errorsion);
    I5out = round(I5out);
    psnr5 = psnr(uint8(round(I5out(:,1:end-offR,:))), Id(:,1:end-offR,:));            
    %ssim5 = ssim(uint8(round(I5out(:,1:end-offR,:))), Id(:,1:end-offR,:));      
    %[psnrhv5,~] = psnrhvsm_mex(uint8(round(I5out(:,1:end-offR,:))), Id(:,1:end-offR,:), 8);
    %PSNRS(v-1, 5) = psnr5;    
    %SSIMS(v-1, 5) = ssim5;    
    %PSNRHVS(v-1, 5) = psnrhv5;    
    figure(a); imshow(uint8(I5out(:,1:end-offR,:))); title(['De With, ',num2str(v-1),'; PSNR=', num2str(psnr5)]); drawnow;
    imwrite(uint8(I5out(:,1:end-offR,:)), ['saves\renders\',dataset,'_',num2str(v),'_dewith.png'], 'png');   
    save(['saves\mat\',dataset,'_',num2str(v),'_dewith.mat'], 'I5out');
end

