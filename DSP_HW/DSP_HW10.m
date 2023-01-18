clear all
%%
im = imread("dunes.jpg"); 
size(im)
figure(1);
imagesc(im)
axis equal
axis tight
colormap gray
% plot a cropped version of the map


imcrop = flipud(im(100:650,650:1200,:));
imBW = rgb2gray(im);
% imBW(imBW<100) = 0;
% imBW(imBW>0) = 1;
imBW = double(im)-mean(double(im(:)));

% assume spacing is 2 m per pixel (close enough!)
kmpp = 2;
easting = (1:size(imBW,2))*kmpp;
northing = (1:size(imBW,1))*kmpp;
imagesc(imBW)
set(gca,'ydir','normal')
axis equal
axis tight
xlabel('easting (m)')
ylabel('northing (m)')
colorbar
figure(2); clf
subplot(2,1,1)
imBW = rgb2gray(im);
% imBW(imBW<1200) = 0;
% imBW(imBW>0) = 1;
imBW = double(imBW)-mean(double(imBW(:)));
easting = (1:size(imBW,2))*kmpp;
northing = (1:size(imBW,1))*kmpp;
imagesc(easting,northing,imBW)
colorbar
xlabel('easting (km)')
ylabel('northing (km)')
set(gca,'ydir','normal')
axis equal
axis tight
colormap gray
% demonstrate how to take an arbitrary profile
x = [0 1200]; % endpoints
y = [0 800]; % endpoints
n = 101; % number of points

xs = linspace(x(1),x(2),n);
ys = linspace(y(1),y(2),n);
hold on
plot(xs,ys,'r.')
% draw a profile of magnetic anomalies
subplot(2,1,2)
mag_prof = interp2(easting,northing,imBW,xs,ys);
dist_prof = sqrt((xs-x(1)).^2 + (ys-y(1)).^2);
plot(dist_prof,mag_prof,'k')
title('magnetic anomalies along cross-section')
xlabel('distance (km)')
axis tight
ylim([-.6 .6])


%% plot 2D amplitude spectrum
figure(3); clf
fftBW = fft2(imBW);
fftshiftBW = (fftshift(fftBW));
ppk = 1/kmpp; % sample rate in pixel per kilometer
subplot(1,3,1) % fft2 amplitude spectrum
fx_axis = linspace(0,ppk,size(fftBW,2)+1);
fx_axis(end) = [];
fy_axis = linspace(0,ppk,size(fftBW,1)+1);
fy_axis(end) = [];
imagesc(fx_axis,fy_axis,abs(fftBW))
xlabel('cycles per mm')
ylabel('cycles per mm')
set(gca,'ydir','normal')
axis equal
axis tight
subplot(1,3,2) % frequency shifted amplitude spectrum
fx_axis_shift = linspace(-ppk/2,ppk/2,size(fftBW,2)+1);
fx_axis_shift(end) = [];
fy_axis_shift = linspace(-ppk/2,ppk/2,size(fftBW,1)+1);
fy_axis_shift(end) = [];
imagesc(fx_axis_shift,fy_axis_shift,abs(fftshiftBW));
xlabel('cycles per m')
ylabel('cycles per m')
set(gca,'ydir','normal')
axis equal
axis tight
subplot(1,3,3) % low-pass filtered amplitude spectrum
alpha = 2; % gaussfilter dimension
smooth_spectrum = imgaussfilt(abs(fftshiftBW),alpha);

imagesc(fx_axis_shift,fy_axis_shift,smooth_spectrum)
xlabel('cycles per m')
ylabel('cycles per m')
set(gca,'ydir','normal')
axis equal
axis tight