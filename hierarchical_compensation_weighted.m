function [Iout, Zout] = hierarchical_compensation_weighted(Iv, Zv, levels)

minZ = min(Zv(:));
maxZ = max(Zv(:));
Iv(repmat(isnan(Zv), [1 1 3])) = nan;
Zv(isnan(Iv(:,:,1))) = nan;
%W{1} = single(~isnan(Zv));
Wx = single((Zv-minZ)/(maxZ-minZ));
Wx(isnan(Wx)) = 0;

I{1} = single(Iv).*repmat(Wx, [1 1 3]);
I{1}(isnan(I{1})) = 0;

Z{1} = single(Zv).*Wx;
Z{1}(isnan(Z{1})) = 0;

W{1} = Wx;

for i=2:levels
    Z{i} = imresize(Z{i-1}, 0.5, 'bilinear');
    I{i} = imresize(I{i-1}, 0.5, 'bilinear');
    W{i} = imresize(W{i-1}, 0.5, 'bilinear');
end

Zout = Z{1};
Iout = I{1};
Wout = W{1};

for i=2:levels
    Z{i} = imresize(Z{i}, size(Zv), 'bilinear');
    I{i} = imresize(I{i}, size(Zv), 'bilinear');
    W{i} = imresize(W{i}, size(Zv), 'bilinear');
    
    Zout = Zout + Z{i}/4^(i-1);
    Iout = Iout + I{i}/4^(i-1);
    Wout = Wout + W{i}/4^(i-1);
end

Zout = single(Zout./Wout);
Zout(~isnan(Zv)) = Zv(~isnan(Zv));
Iout = single(Iout./repmat(Wout, [1 1 3]));
Iout(~isnan(Iv)) = Iv(~isnan(Iv));

