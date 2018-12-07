I = single(imread('view1.png'));
D = single(imread('disp1_noholes.png'));
D(D<60) = nan;
[h w] = size(D);
HW = h * w;
figure; imshow(D, []); colormap(pink);
%%
A = D(1430:end, :);
A(A<130) = nan;
D(1430:end, :) = A;
figure; imshow(D, []); colormap(pink); title('With Holes')
D2 = D;
D2(isnan(D)) = 0;
%%
D3 = recursive_bilateral(D2, I)./recursive_bilateral(single(~isnan(D)), I);
D3(~isnan(D)) = D(~isnan(D));
figure; imshow(D3, []); colormap(pink); title('Compensated')
imwrite(uint8(round(D3)), 'disp1_noholes.png', 'png');
%%
L = get_laplacian(I);
%%
S = spdiags(double(~isnan(D(:))), 0, HW, HW);
%%
lambda = 1e-4;
Y = D(:);
Y(isnan(Y)) = 0;
X = (S'*S + lambda*L)\S'*Y;
figure; im);show(reshape(X, [h w]), []);