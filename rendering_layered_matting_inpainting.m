function [Iv] = rendering_layered_matting_inpainting(Cv, I, Z, C, errosion_radius, avg_alpha, hf_alpha)

r = errosion_radius; % errosion radius
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

%figure; imshow(imerode(T2==2, strel('disk', r)) - imerode(T2==0, strel('disk', r)), [-1 1]); title('Trimap 2');


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
    %figure; imshow(uint8(B2i)); title('Recursive-Compensated BG');

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