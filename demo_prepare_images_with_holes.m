close all
clear
clc

dataset = 'Aloe'; 
%dataset = 'Art'; 
%dataset = 'Cones'; 
%dataset = 'Reindeer'; 
%dataset = 'Moebius'; 
%dataset = 'Teddy';



I1 = single(imread(['datasets\',dataset,'\view1.png']));
D1 = single(imread(['datasets\',dataset,'\disp1_noholes.png']));

[h, w] = size(D1);
%HW = h * w;
maxdisp = ceil(max(D1(:)));
mindisp = 0;
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

%
v=2; % desired view index
offR = 50*(v-1);
t = [-40*(v-1) 0 0]';
%t = [-50 -50 1000]';
Cv = single(K*[eye(3), t]);    
radius = 0;
%
[Iv, Zv, ~] = free_rendering(Cv, I1, Z1, C1);
if radius > 0
    Valid = imerode(~isnan(Zv),strel('disk',radius));
    Iv(repmat(~Valid, [1 1 3])) = nan;
    Zv(~Valid) = nan;
end
Iholes = uint8(Iv(:,1:end-offR,:));
figure; imshow(Iholes); title('Rendered with no filling');
imwrite(Iholes, ['saves/',dataset,'_',num2str(v),'_holes.jpg'], 'jpg', 'Quality', 80);
%imwrite(Iholes, ['saves/',dataset,'_',num2str(v),'_holes.jpg'], 'jpg', 'Quality', 80);
%imwrite(Iholes, ['saves/',dataset,'_',num2str(v),'_right_holes.jpg'], 'jpg', 'Quality', 80);
%imwrite(Iholes, ['saves/',dataset,'_',num2str(v),'_holes_4.jpg'], 'jpg', 'Quality', 80);


