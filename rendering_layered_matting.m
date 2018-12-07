function [Iv, Ar, Ifr, Ibr, If, Ib, Zf, Zb, alpha] = rendering_layered_matting(Cv, I1, Z1, D1, C1, errosion, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda)
%errosion; % errosion radius
%avg_alpha = 0.96;
%avg_alpha2 = 0.98;
alpha_thr = 0.01;
t_thr = 3;
[h, w] = size(Z1);
HW = h * w;
D_hat = recursive_gaussian(D1, avg_alpha)./recursive_gaussian(ones([h w], 'single'), avg_alpha);
FgA = D1 > (D_hat + t_thr);

% some overlap between Bg and Fg is required
Both = imopen((D1 > D_hat) & ~(D1 > D_hat + 3), strel('disk', t_thr));

FgDepthMask = FgA | Both;
BgDepthMask = ~FgA | Both;
%BgDepthMask = imerode(~(D1 > D_hat), strel('disk', 1)) | Both;
%BgDepthMask = ~(D1 > D_hat) | Both;


%FgB = D1 > (D_hat + t_thr);
%BgB = imerode(~FgB, strel('disk', errosion));

%figure; imshow(FgA); title('Foreground A');
%figure; imshow(BgB); title('Background B');

%[Dx, Dy] = gradient(D1);
%W = ~(D1 > D_hat);
%W = single(imerode(W, strel('disk', errosion)));
%figure; imshow(W); title('Weight');

T2 = FgA + 1 - (~FgA);

%figure; imshow(imerode(T2==2, strel('disk', r)) - imerode(T2==0, strel('disk', r)), [-1 1]); title('Trimap 2');


    Ti = uint8(1 + imerode(T2==2, strel('disk', errosion)) - imerode(T2==0, strel('disk', errosion)));
    %figure; imshow(Ti, [0 2]); title('Trimap');
    clear mex
    %tic
    %alpha = learning_matting(I1, Ti, 1, 10);
    %toc

    AL = get_learning_laplacian(uint8(I1), Ti, matting_radius, matting_sigma, matting_epsilon);
    S = double(Ti(:)~=1);
    S = spdiags(S(:), 0, HW, HW);
    alpha = (matting_lambda*AL + S'*S) \ (S'*(double(Ti(:))));
    alpha = single(reshape(alpha, [h w]))/2;
    %alpha = single(reshape(alpha, [h w]) + 1)/2;
    alpha(alpha < 0) = 0;
    alpha(alpha > 1) = 1;
    %figure; imshow(alpha, [0 1]); title('Reconstructed Alpha');
    [Fi, Bi] = solveFB(double(I1)/255, double(alpha));
    Fi = single(Fi*255);
    Bi = single(Bi*255);    
    Fi(repmat(Both, [1 1 3])) = I1(repmat(Both, [1 1 3]));
    Bi(repmat(Both, [1 1 3])) = I1(repmat(Both, [1 1 3]));
    
    %figure; imshow(uint8(Fi)); title('Recovered Foreground');
    %figure; imshow(uint8(Bi)); title('Recovered Background');
    
    FgMask = (alpha >= 0.2) | Both;
    
    %FgMask = (alpha >= (1-alpha_thr)) | Both;
    %BgMask = (alpha <= alpha_thr) | Both;
    % simple compensation
    %If = I1;
    If = Fi;
    %F2i(repmat(alpha < 0.5, [1 1 3])) = nan;    
    If(repmat(~FgMask, [1 1 3])) = nan;    
    If = recursive_compensation_simple(If, hf_alpha);
    %figure; imshow(uint8(If)); title('Compensated Foreground');
    
    %B2i = Bi;
    %B2i(repmat(alpha > 0, [1 1 3])) = nan;
    
    %BgMask = (alpha <= 0.9) | Both;
    %Ib = Bi;
    BgMask = (alpha <= 0.05) | Both;
    Ib = I1;
    Ib(repmat(~BgMask, [1 1 3])) = nan;
    Ib = recursive_compensation_simple(Ib, hf_alpha);
    %figure; imshow(uint8(Ib)); title('Compensated Background');

    
    % compensation of depth maps
    Zf = Z1;
    Zf(~FgDepthMask) = 0;        
    Zf = recursive_gaussian(Zf, hf_alpha)./recursive_gaussian(single(FgDepthMask), hf_alpha);
    Zf(FgDepthMask) = Z1(FgDepthMask);    
    %Df = recursive_compensation_simple(Df, hf_alpha);
    %Zf = f*baseline./(Df);
    clear Df;
    
    Zb = Z1;
    Zb(~BgDepthMask) = 0;        
    Zb = recursive_gaussian(Zb, hf_alpha)./recursive_gaussian(single(BgDepthMask), hf_alpha);  
    %Db = recursive_compensation_simple(Db, hf_alpha);        
    Zb(BgDepthMask) = Z1(BgDepthMask);    
    %Zb = f*baseline./(Db);
    clear Db;
        
    
    [Ifr, ~, ~] = free_rendering(Cv, If, Zf, C1);
    [Ar, ~, ~] = free_rendering(Cv, alpha, Zf, C1);    
    [Ibr, Zbr, ~] = free_rendering(Cv, Ib, Zb, C1); % blurry
    %figure; imshow(isnan(B2v(:,:,1))); title('Bg Nans')    
    %
    %[Ibr] = recursive_compensation(Ibr, Zbr, 0.6);
    [Ibr] = recursive_compensation_simple(Ibr, 0.6);
    
    %figure; imshow(uint8(Ibr)); title('Compensated Bg')    
    Ifr(isnan(Ifr)) = 0;
    Ar(isnan(Ar)) = 0;
    
    Iv = Ifr.*repmat(Ar, [1 1 3]) + Ibr.*repmat((1-Ar), [1 1 3]);      