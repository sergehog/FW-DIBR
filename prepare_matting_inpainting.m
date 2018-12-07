function [alpha, If, Ib, Df, Db] = prepare_matting_inpainting(I1, D1, errosion, avg_alpha, hf_alpha, matting_radius, matting_sigma, matting_epsilon, matting_lambda)
alpha_thr = 0.01;
t_thr = 3;
[h, w] = size(D1);
HW = h * w;
D_hat = recursive_gaussian(D1, avg_alpha)./recursive_gaussian(ones([h w], 'single'), avg_alpha);
FgA = D1 > (D_hat + t_thr);

% some overlap between Bg and Fg is required
%Both = imopen( (D1 > D_hat)  & ~(D1 > D_hat + 3), strel('disk', t_thr));
Both = (D1 > D_hat)  & ~(D1 > D_hat + 3);

FgDepthMask = FgA | Both;
%BgDepthMask = imerode(~(D1 > D_hat), strel('disk', 1)) | Both;
BgDepthMask = ~(D1 > D_hat) | Both;

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
    %If = recursive_compensation_simple(If, hf_alpha);
    %figure; imshow(uint8(If)); title('Compensated Foreground');
    
    %B2i = Bi;
    %B2i(repmat(alpha > 0, [1 1 3])) = nan;
    
    %BgMask = (alpha <= 0.9) | Both;
    %Ib = Bi;
    BgMask = (alpha <= 0.0) | Both;
    Ib = I1;
    Ib(repmat(~BgMask, [1 1 3])) = nan;
    %Ib = recursive_compensation_simple(Ib, hf_alpha);
    %figure; imshow(uint8(Ib)); title('Compensated Background');
    
 % compensation of depth maps
    Df = D1;
    Df(~FgDepthMask) = nan;        
    %Df = recursive_gaussian(Df, hf_alpha)./recursive_gaussian(single(FgDepthMask), hf_alpha);
    %Df(FgDepthMask) = D1(FgDepthMask);    
    %Df = recursive_compensation_simple(Df, hf_alpha);
    %Zf = f*baseline./(Df);
    %clear Df;
    
    Db = D1;
    Db(~BgMask) = nan;    
    
    %Db = recursive_gaussian(Db, hf_alpha)./recursive_gaussian(single(BgDepthMask), hf_alpha);  
    %Db = recursive_compensation_simple(Db, hf_alpha);        
    %Db(BgDepthMask) = D1(BgDepthMask);    
    %Zb = f*baseline./(Db);
    %clear Db;    
