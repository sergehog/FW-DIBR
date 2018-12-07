clc
close all
clear 

%dataset = 'Aloe'; x1 = 400; x2 = 1032; y1 = 200; y2 = 900;
dataset = 'Art'; x1 = 400; x2 = 1032; y1 = 200; y2 = 900;

%cd 'saves\renders';
files = dir;
af = figure();
for i=3:size(files)
    if strfind(files(i).name, dataset) > 0
        I = imread(files(i).name);
        I = I(y1:y2, x1:x2, :);
        figure(af);
        imshow(I);
        title(files(i).name);
        drawnow;
        %break;
        imwrite(I, ['..\',dataset,'\',files(i).name], 'jpg');
    end
end
%%
