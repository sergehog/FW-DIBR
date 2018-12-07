clear
close all
clc

%dataset = 'Art'; 
%dataset = 'Reindeer'; 
%dataset = 'Aloe'; 
%dataset = 'Moebius'; 
%dataset = 'Cones'; 
dataset = 'Teddy';

%I1 = single(imread(['datasets\',dataset,'\view1.png']));
%D1 = single(imread(['datasets\',dataset,'\disp1_noholes.png']));
%[h, w] = size(D1);
%mindisp = floor(min(D1(:)));
%maxdisp = ceil(max(D1(:)));
%addpath(genpath('../matting'));

for i=1:6
    Is{i} = imread(['datasets\',dataset,'\view',num2str(i),'.png']);
end

algos = {'hhf', 'jtdi', 'dewith', 'fukushima', 'layered', 'proposed'};
names = {'DA-HHF', 'JTDI', 'DGI-FFV', 'BLUR', 'Layered Simple', 'Layered Inpainting'};
%
clc
PSNRS = zeros([5 numel(algos)]);
PSNRHVS = zeros([5 numel(algos)]);
SSIMS = zeros([5 numel(algos)]);
af = figure();
for v=2:6  
    offR = 50*(v-1);
    Id = Is{v};
    Id = Id(:,1:end-offR,:);
       
    [h w ~] = size(Id);

    if (strcmp(dataset,'Teddy') ~= 0 || strcmp(dataset,'Cones') ~= 0) && v==6        
       Id = Id(:,53:end,:);        
    end

    for a=1:numel(algos)
        algo = algos{a};
        Iout = imread(['saves\renders\',dataset,'_',num2str(v),'_',algo,'.png']);
        Iout = Iout(1:h, 1:w, :);
        
        if (strcmp(dataset,'Teddy') ~= 0 || strcmp(dataset,'Cones') ~= 0) && v==6        
            Iout = Iout(:,53:end,:);
        end
    
        psnrI = psnr(uint8(round(Iout)), Id);            
        ssimI = ssim(uint8(round(Iout)), Id);      
        [psnrhvI,~] = psnrhvsm_mex(uint8(round(Iout)), Id, 8);
        PSNRS(v-1, a) = psnrI;    
        SSIMS(v-1, a) = ssimI;    
        PSNRHVS(v-1, a) = psnrhvI;    
        %figure(af); imshow(uint8(Iout)); 
        %title([algo, '; PSNR=',num2str(psnrI), '; SSIM=', num2str(ssimI), '; PSNRHVS=', num2str(psnrhvI)]);
        %drawnow;
    end
end

%
%PSNRS = PSNRS(:,1:4)

fig1 = figure; 
plot(PSNRS, 'LineWidth', 2); 
legend(names);
title(['Rendering performance for "',dataset,'" dataset']);
ylabel('PSNR, dB');
xlabel('Rendering offset');
%ylim([25 32]);
%ylim([22 33]);
saveas(fig1, [ 'figures/',dataset,'_psnr.pdf'], 'pdf');
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
