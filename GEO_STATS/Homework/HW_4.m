%% HW  4 
close all; clear all
%%  Problem 1 
D = load('DevonBdot.txt'); 
dist = D(:,1); 
bdot = D(:,2);
figure(1);clf;
clear D2 D4
plot(dist, bdot); 
D1 = D(1:100,:);
D2 = D(101:200,:); 
D3 = D(201:300,:); 
D4 = D(301:400,:); 
D1_mu = mean(D1);
D2_mu = mean(D2);
D3_mu = mean(D3); 
D4_mu = mean(D4);
D1_sig = std(D1); 
D2_sig = std(D2); 
D3_sig = std(D3); 
D4_sig = std(D4); 
%calculate vairance for rates
D1_var = var(D1);
D2_var = var(D2);
D3_var = var(D3); 
D4_var = var(D4);

%% plot relative density plots 
figure(1); clf 
subplot(2,2,1);
x0 = 15; %define number of bins ;
[distr, xbins]= hist(D1(:,2), x0);
dx = xbins(2) -xbins(1); %bin width by subtracting bin centers (xbins)
distr_area = sum(distr*dx);
norm =distr/distr_area; 
bar(xbins, norm); hold on  %plot the relative  density histrogram problem

% kernal density histrogram 
h = 1; 
xm = D1_mu(2)-10*D1_sig(2):D1_sig(2)/10:D1_mu(2)+10*D1_sig(2);
for n=1:length(xm)

    dist = D1(:,2)-xm(n); %distance from x0 to all other data 
    Ix = find(abs(dist)<h) % finding data within window h
    w =15/16*(1-dist(Ix)/h).^2
    f1(n) = sum(w); %store the estiamte for each position in our array  
end 
dw = D1_sig(2)/10; %width of bins 
f1 = 1/sum(f1.*dw)*f1 ;%normalizes the pdf 

plot(xm, f1, 'r', 'linewidth', 2)
title( '0m-200m')
% section2 
subplot(2,2,2);
x0 = 15; %define number of bins ;
[distr, xbins]= hist(D2(:,2), x0);
dx = xbins(2) -xbins(1); %bin width by subtracting bin centers (xbins)
distr_area = sum(distr*dx);
norm =distr/distr_area; 
bar(xbins, norm); hold on  %plot the relative  density histrogram problem

% kernal density histrogram 
h = 1; 
xm = D2_mu(2)-10*D2_sig(2):D2_sig(2)/10:D2_mu(2)+10*D2_sig(2);
for n=1:length(xm)

    dist = D2(:,2)-xm(n); %distance from x0 to all other data 
    Ix = find(abs(dist)<h) % finding data within window h
    w =15/16*(1-dist(Ix)/h).^2
    f1(n) = sum(w); %store the estiamte for each position in our array  
end 
dw = D2_sig(2)/10; %width of bins 
f1 = 1/sum(f1.*dw)*f1 ;%normalizes the pdf 

plot(xm, f1, 'r', 'linewidth', 2)
title( '200m-400m')
%Section 3
subplot(2,2,3);
x0 = 15; %define number of bins ;
[distr, xbins]= hist(D3(:,2), x0);
dx = xbins(2) -xbins(1); %bin width by subtracting bin centers (xbins)
distr_area = sum(distr*dx);
norm =distr/distr_area; 
bar(xbins, norm); hold on  %plot the relative  density histrogram problem

% kernal density histrogram 
h = 1; 
xm = D3_mu(2)-10*D3_sig(2):D3_sig(2)/10:D3_mu(2)+10*D3_sig(2);
for n=1:length(xm)

    dist = D3(:,2)-xm(n); %distance from x0 to all other data 
    Ix = find(abs(dist)<h) % finding data within window h
    w =15/16*(1-dist(Ix)/h).^2
    f1(n) = sum(w); %store the estiamte for each position in our array  
end 
dw = D3_sig(2)/10; %width of bins 
f1 = 1/sum(f1.*dw)*f1 ;%normalizes the pdf 

plot(xm, f1, 'r', 'linewidth', 2)
title( '400m-600m')
%Section4 
subplot(2,2,4);
x0 = 15; %define number of bins ;
[distr, xbins]= hist(D4(:,2), x0);
dx = xbins(2) -xbins(1); %bin width by subtracting bin centers (xbins)
distr_area = sum(distr*dx);
norm =distr/distr_area; 
bar(xbins, norm); hold on  %plot the relative  density histrogram problem

% kernal density histrogram 
h = 1; 
xm = D4_mu(2)-10*D4_sig(2):D4_sig(2)/10:D4_mu(2)+10*D4_sig(2);
for n=1:length(xm)

    dist = D4(:,2)-xm(n); %distance from x0 to all other data 
    Ix = find(abs(dist)<h) % finding data within window h
    w =15/16*(1-dist(Ix)/h).^2
    f1(n) = sum(w); %store the estiamte for each position in our array  
end 
dw = D4_sig(2)/10; %width of bins 
f1 = 1/sum(f1.*dw)*f1 ;%normalizes the pdf 

plot(xm, f1, 'r', 'linewidth', 2)
title( '600m-800m')
%% Problem 2 
% the answers to the questions are that the mean for the acumulation rates
% is fairly constant but the the variance is not constant 
%% Problem 3 ? calculate semivariance 
clear dist ; 
dist = D(:,1); 
bdot = D(:,2);
[h,V,c] = semivariogram(dist, bdot);
figure(2); clf;
plot(h,V,'o');hold on 
title('Semivariance and Covariance plots')%plots semivariance 
plot(h,c, 'o'); %plots covariance 
% plot autocorrelation 

figure(3); clf 
autocorr(bdot, 'NumLags', 399); 


%% problem 4 
 [h,V,npairs] = semivariogram_mc(dist, bdot, 10); 
 clear figure(3);
 figure(3); hold on 
 plot(h, V, 'k'); hold on
 [h,V,npairs] = semivariogram_mc(dist, bdot, 50); 
  plot(h, V, 'r'); hold on
 [h,V,npairs] = semivariogram_mc(dist, bdot, 100); 
  plot(h, V, 'b'); hold on
 legend('npair =10', 'npair=50', 'npair = 100')
   
 figure(4); hold on; clf 
    plot(h, npairs); hold on
%% Problem 5 Plot semivariance for 4 seperate sections 

[h,V,c] = semivariogram(D1(:,1), D1(:,2));
figure(5); clf;
plot(h,V,'o', 'LineWidth', 2);hold on 
title('Semivariance plots')%plots semivariance 

[h,V,c] = semivariogram(D2(:,1), D2(:,2));
plot(h,V,'o', 'LineWidth', 2);hold on 
title('Semivariance plots')%plots semivariance 

[h,V,c] = semivariogram(D3(:,1), D3(:,2));
plot(h,V,'o', 'LineWidth', 2);hold on 
title('Semivariance plots')%plots semivariance 

[h,V,c] = semivariogram(D4(:,1), D4(:,2));plot(h,V,'o','LineWidth', 2);hold on 
title('Semivariance plots')%plots semivariance 
legend('Section 1', 'Section 2', 'Section 3', 'Section 4')
%% Answer to Question 
% Yes, semivariance relies on the lag. The sill is the range of the
% measured
% semivariance value which depeneds on the square of the standard
% deviations and the standard deviation is the range of error you can
% expect for a given value. The standard deviation is used to determine the
% semivariance therfore making the semivarience dependenent on that range
% of error or lag. 
%% Problem 6 finding best vaarigram parameters 

a = 1:60 ; %range of a values 
c = 30:60 ; %range of sill values 
for p =1:length(a) %loop through values of bdot
    for q=1:length(c) %loop through values of A 
     rmse(p,q)= model_variogram(h,V,c(q),a(p),'L');
    end 
end 
figure(3); clf 
imagesc(c,a,rmse);
colorbar
% V = model_variogram(h,c,a,n,type)
%% Problem 7 repeat above using fminsearch 
fh = @(p) model_variogram(h,V,p(1), p(2), 'L');
[pbest, fval] = fminsearch(fh,[30,60]);
figure(3); hold on
plot(pbest(1), pbest(2), 'wo', 'markersize', 8, 'linewidth', 2) 
%% Probelm 8 
fh = @(p) model_variogram(h,V,p(1), p(2), 'S');
[pbest, fval] = fminsearch(fh,[30,60]);
figure(3); hold on
plot(pbest(1), pbest(2), 'wo', 'markersize', 8, 'linewidth', 2) 