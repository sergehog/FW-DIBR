function M = local_otsu_thresholding(S, radius)
S = S./max(S(:));
[h, w] = size(S);
M = zeros([h w], 'uint8');
for x=1:w
    for y=1:h
        Block = S(max(1, y-radius):min(h, y+radius), max(1, x-radius):min(w, x+radius));        
        %M(max(1, y-radius):min(h, y+radius), max(1, x-radius):min(w, x+radius)) = M(max(1, y-radius):min(h, y+radius), max(1, x-radius):min(w, x+radius)) + uint8(Block > graythresh(Block));
        M(y,x) = S(y,x) > graythresh(Block);
    end
end