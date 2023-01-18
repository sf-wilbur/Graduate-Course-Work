close all; clear all

%% 1: Enter and run the script below. Change the values of sigma and
%subplots 3 and 4 you will have to design a complementary low-pass
%filter (using sigma=5) and plot it and its frequency response.
%Recall that the way to create a high pass filter is to add a central
%delta function to the negative values of the low pass FIR filter.

figure(1); subplot(2,2,1) % low-pass gaussain spatial filter
sigma = 5;
m = 17; n = 17;
flow = fspecial('gaussian',[n m],[sigma]);
mesh(0:m-1,0:n-1,flow);
title('gaussian 2-D low pass filter')
subplot(2,2,2) % and its frequency response
flow_fft = abs(fft2(flow));
mesh(0:m-1,0:n-1,flow_fft);
title('gaussian 2-D filter amplitude response')


subplot(2,2,3) % create high pass filter from the gaussian

flow_high = -flow;
flow_high(9,9) = flow_high(9,9)+1;
mesh(0:m-1,0:n-1,flow_high)
title('high pass filter gaussian')


subplot(2,2,4) % derive the amplitude response from above
title('high pass filter amplitude response')
fhigh_noshift = abs(fft2(flow_high));
fhigh_fft = fftshift(fhigh_noshift); 
mesh(0:m-1,0:n-1,fhigh_fft)
title('high pass filter amplitiude response')
%% Problem 2 :
load('Idaho_grav_grid.mat')
load('grav_data.mat')

figure(2) 
imagesc(lons,lats,grav);hold on
plot(Idaho.Lon, Idaho.Lat,'k', 'LineWidth',2)
set(gca, 'ydir', 'normal')
daspect([1 1 1])
c = colorbar;

%% Problem 3 
amp = abs(fft(grav));
gspec = log10(amp);
f1 = 0:length(grav(:,1)-1);
f2 = 0:length(grav(1,:)-1);
subplot(2,2,1) % 2D fft
figure(3)
imagesc(f1,f2, gspec)
title('FFT Gravity Data ')
xlabel('East West')
ylabel('North to South')
xlim([0 125]) % this is the nyquist (east to west)
ylim([0 75]) %this is nyquist for lat (north to south)
colorbar


%%
lp = conv2(flow, gspec);
hp = conv2(flow_high, gspec);

figure(4); clf 
subplot(2,1,1)
imagesc(f1,f2,lp); % low pass
title('Low Pass Filtered FFT')
xlabel('Cycles per km by Latitude')
ylabel('Cycles per km by Longitude')

xlim([1 125])
ylim([0 75])


subplot(2,1,2)

imagesc(f1,f2,hp); %high pass
title('High Pass Filtered FFT')
xlabel('Cycles per km by Latitude')
ylabel('Cycles per km by Longitude')
xlim([1 125])
ylim([0 75])


subtitle('Low pass and high pass filter of Idaho Gravity Map')

