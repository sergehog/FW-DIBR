function [psnr_value] = PSNR(I1, I2)
%I1 = single(I1);
%I2 = single(I2);
%I1(isnan(I1)) = 0;
%I2(isnan(I2)) = 0;
rmse_value = sqrt(mean((I1(:)-I2(:)).^2));
psnr_value = 20*log10(255/rmse_value);