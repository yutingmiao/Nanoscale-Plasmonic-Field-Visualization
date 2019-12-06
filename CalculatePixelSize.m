%% Read the .tif file
clear all; close all;
f = '6';    
fname = [f, '.tif'];    % filename
im = double(imread(fname));
figure;imshow(im,[]);

% Gaussian filter to smooth image (optional)
% d = zeros(11,11); d(1,1)=-1; d(1,2)=1;    
% H = fspecial('gaussian', [13,13], 1);
% ims = imfilter(im,H);    
% dim = abs(conv2(ims,d));

%% Set parameters
x = [113, 427];    % crop out area of the structure with known dimensions (x and y are converted)
y = [124, 434];
figure(1);imshow(im(x(1):x(2), y(1):y(2)),[]);

thresX = 45000; thresY = 45000;    % set the threshold for the peak selection
disX = 200; disY = 200;

%% Peak selection
locx_all = []; locy_all = []; px_all = []; py_all = [];
for i0 = x(1):x(2)
    im_x = im(i0, y(1):y(2));
    [px, locx] = findpeaks(im_x,'MinPeakHeight',thresX,'MinPeakDistance',disX);
    warning('off','signal:findpeaks:largeMinPeakHeight')
    figure(2); hold on; plot(im_x);plot(locx,px,'o');
    locx_all = [locx_all, locx]; px_all = [px_all, px];
end; title('hori')

for j0 = y(1):y(2)
    im_y = im(x(1):x(2), j0);
    [py, locy] = findpeaks(im_y,'MinPeakHeight',thresY,'MinPeakDistance',disY);
    warning('off','signal:findpeaks:largeMinPeakHeight')
    figure(3); hold on; plot(im_y);plot(locy,py,'o');
    locy_all = [locy_all, locy.']; py_all = [py_all, py.'];
end; title('vert')

%% Plot figures and calculate pixel size
figure(4); hold on; plot(locx_all,'x'); title('Hori Edge');
figure(5); hold on; plot(locy_all,'x'); title('Vert Edge');

locx_all = locx_all(24:475);    % select only data points with sharp edges 
figure(6); hold on; plot(locx_all,'x'); title('Hori Edge');
left = locx_all < 50;    % select data points only for left edge
left=locx_all(left);
right = locx_all > 250;
right=locx_all(right);

locy_all = locy_all(211:484);
figure(7); hold on; plot(locy_all,'x'); title('Vert Edge');
bottom = locy_all < 50;
bottom=locy_all(bottom);
top = locy_all > 200;
top=locy_all(top);

pnum_tb = mean(top) - mean(bottom);    % calculate the pixel number
pnum_lr = mean(right) - mean(left);

ssize_tb = 20000; ssize_lr = 20000;
psize = (ssize_tb / pnum_tb + ssize_lr / pnum_lr) / 2;
% 82-83 nm/pixel