% Geostats: 09/16/2021
clear all;
D = load('elevations.txt');

D = D(:); % creates one big long column 
Ix = find(isfinite(D)) %isfinits is putting zeros in a new vector where there are non-finite values 
                        % and find is finding the values where there are ones and not zeros 

D = D(Ix)

     % and find is finding the values where there are ones and not zeros 
%% PROBLEM 1

D_min = min(D);
D_max = max(D);
D_std = std(D);
D_mean= mean(D);
nbins = 100;
x0 = (D_mean-100*D_std):D_std/100:(D_mean+100*D_std); %sort the values from least to greatest by 10 off the average

%% PROBLEM 2
figure(1)
 
[distr, xbins]= hist(D, x0);
dx = D_std/10; %bin width
distr = distr / sum(distr*dx);
bar(xbins, distr) %plot the relative  density histrogram problem

%% PROBLEM 3
RD = randsample(D,10,true) % draws 10 random numbers from D, the true removes the number and stores it and picks remaining values from new list of data
RD_mean = mean(RD)
RD_min= min(RD)
RD_max= max(RD)

%% PROBLEM 4

for i=1:1000;

   RD = randsample(D,10,true);
   RD_mean(i)= mean(RD);
   RD_min(i) = min(RD);
   RD_max(i) = max(RD);
   RD_std(i) = std(RD);
  
end
%% PROBLEM 5
figure(2)
[distr_mean, xbins_mean]= hist(RD_mean, nbins);
dx_mean = xbins_mean(2)-xbins_mean(1);
distr_mean = distr_mean / (sum(distr_mean.*dx_mean));
bar(xbins_mean, distr_mean) %plot the relative  density histrogram problem

figure(3)
[distr_std, xbins_std]= hist(RD_std, nbins);
dx_std = xbins_std(2)-xbins_std(1);
distr_std = distr_std / sum(distr_std);
bar(xbins_std, distr_std) %plot the relative  density histrogram problem

figure(4)
[distr_min, xbins_min]= hist(RD_min, nbins);
dx_min = xbins_min(2)-xbins_min(1);
distr_min = distr_min / sum(distr_min.*dx_min);
bar(xbins_min, distr_min) %plot the relative  density histrogram problem


figure(5)
[distr_max, xbins_max]= hist(RD_max, nbins);
dx_max = xbins_max(2)-xbins_max(1);
distr_max = distr_max / sum(distr_max.*dx_max);
bar(xbins_max, distr_max) %plot the relative  density histrogram problem


%% PROBLEM 6
ND = load('elevations.txt');
ND = ND(:);
zx = std(D)/10;
z0 = (mean(D)-std(D)*10):zx:(mean(D)+std(D)*10)
[N, xbins]= hist(D,z0);
RDH = N/sum(N*zx);
figure(6);
bar(z0,RDH)
h = 10; %window size 
clear f;
for n=1:length(z0)
    dist = (D-z0(n)); %distance from x0 to all other data values
    Idx = find(abs(dist)<h); %finding all data points within h from x0
    w =15/16*(1-(dist(Idx)/h).^2).^2; %weights of all points
    f(n) =sum(w); %sum the weights 

end 

zx = std(D)/10;
f= 1/sum(f*zx)*f;%normalized pdf 
gauss = mypdf(z0,D_mean, D_std);


figure(6); hold on 
 plot(z0, f, 'g', 'linewidth', 2)
plot(z0, gauss, 'r')

%% PROBLEM 7
%% for ther sample mean
figure(7); hold on
bar(xbins_mean, distr_mean)%For mean distribution values 
mu_0 = std(RD_mean)/10;
RD_0 = (mean(RD_mean)-std(RD_mean)*10):mu_0:(mean(RD_mean)+std(RD_mean)*10)
RDmu = mean(RD_mean)
RDsig = std(RD_mean)
mu_g = mypdf(RD_0, RDmu, RDsig)

plot(RD_0, mu_g, 'r', 'linewidth', 2)
%% for the sample standard deviation
figure(8); hold on
bar(xbins_std, distr_std)
std_0 =std(RD_std)/10;
RDstd_0 = (mean(RD_std)-std(RD_std)*10):std_0:(mean(RD_std)+std(RD_std)*10)
STDmu = mean(RD_std)
STDsig = std(RD_std)
STD_g = mypdf(RDstd_0, STDmu, STDsig)
plot(RDstd_0, STD_g, 'r', 'linewidth', 2)


%% for the sample max 
clear f;
figure(9); hold on
bar(xbins_max, distr_max);
max0 =std(RD_max)/10;
Max_0 = (mean(RD_max)-std(RD_max)*10):max0:(mean(RD_max)+std(RD_max)*10);
Maxmu = mean(RD_max);
Maxsig = std(RD_max);
max_g = mypdf(Max_0, Maxmu, Maxsig);
plot(Max_0, max_g, 'r', 'linewidth', 2)

%% for sample min

clear f;
figure(10); hold on
bar(xbins_min, distr_min);
max0 =std(RD_min)/10;
Max_0 = (mean(RD_min)-std(RD_min)*10):max0:(mean(RD_min)+std(RD_min)*10);
Maxmu = mean(RD_min);
Maxsig = std(RD_min);
max_g = mypdf(Max_0, Maxmu, Maxsig);
plot(Max_0, max_g, 'r', 'linewidth', 2)
%% question 8 

IDmu = find(RD_mean<=D_mean) 
perc_mu = length(IDmu)/length(RD_mean) 
perc_mu = perc_mu*100

fprintf( '%f chance of finding a value less than the true mean', perc_mu) 

%% Question 9

tru_max = find(RD_max>(0.99*D_max) & RD_max<(1.01*D_max)) 
tru_min = find(RD_min>(0.99*D_min) & RD_min<(1.01*D_min)) 


%range = tru_max(min(tru_max):max(tru_max)):tru_min(min(tru_min):max(tru_min):2)
perc_max = length(tru_max)/length(RD_max)
perc_min  = length(tru_min)/length(RD_min)
%Should this be one probabilty 

%% Question 10 
Elev_Range = find(RD_mean >(D_mean - RDsig) & RD_mean < (D_mean + RDsig))
Prob = RD_mean(Elev_Range)
Prob_min = min(Prob) 
Prob_max = max(Prob)

%% question 11
clear D;
D = load('elevations.txt');

D = D(:); % creates one big long column 
Ix = find(isfinite(D)) %isfinits is putting zeros in a new vector where there are non-finite values 
                        % and find is finding the values where there are ones and not zeros 

D = D(Ix)

samplesize = [10, 50, 100, 500];
for q= 1:length(samplesize)
    for n = 1:1000
      D2 = randsample(D, samplesize(q), true);
      Dmean(n,q) = mean(D2); % mean of D2
      Dmin(n,q) = min(D2) ;
      Dstd(n,q) = std(D2);
      Dmax(n,q) = max(D2) ;
    end 
end
figure(14)
subplot(2,2,1)
boxplot(Dmean)
title('D Mean')
subplot(2,2,2)
boxplot(Dstd)
title('D Standard Deviation')
subplot(2,2,3)
boxplot(Dmax)
title('D Maximum')
subplot(2,2,4)
boxplot(Dmin)
title('D Minimum')
%% 
%question 12 Uniform Sampling 
clear D,Ix, D2, x0, D_std; 
figure(11)
D = load('elevations.txt')
nx = 22
nr = 111
nc = 110
Ix=1:nx:nc
Iy = 1:nx:nr
% this time we leave D as a matrix (110x111)
%plot the realtive density histogram 
D_idx = zeros(size(D));
%index D with a vector that describes the rows and clolumns we want 
Ix = find(isfinite(D)) %isfinits is putting zeros in a new vector where there are non-finite values 
                        % and find is finding the values where there are ones and not zeros 

for x = 1:1:22
    for y= 1:1:22
    Iy1 = y:nx:nr;
    Ix1 = x:nx:nc;
    D_idx = D(Ix1, Iy1);
    D_idx = D_idx(isfinite(D_idx));
    D_idx = D_idx(:);
    Didx_mu_200(x,y) = mean(D_idx);
  
    Didx_min_200(x,y) = nanmin(D_idx(:));
    Didx_max_200(x,y) = nanmax(D_idx(:));
    Didx_std_200(x,y) = nanstd(D_idx(:));
    end 
end
% histogram(Didx_mu_200, 30,'Normalization','pdf')

D_mu = mean(Didx_mu_200, 'all');
D_std = std(Didx_mu_200(:));
x0 =(D_mu-100*D_std):D_std/100:(D_mu+100*D_std)
figure
histogram(Didx_std_200, 30,'Normalization','pdf'); hold on 
% gauss_c = mypdf(x0, D_mu, D_std);
% plot(x0, gauss_c, 'r') 



%% Question 13
figure(12)
clear D, D2, x0, D_std, D_mu, Ix1, Iy1, D_idx,D_mu; 

D = load('elevations.txt')
nx = 3
nr = 111
nc = 110
Ix=1:nx:nc
Iy = 1:nx:nr
% this time we leave D as a matrix (110x111)
%plot the realtive density histogram 
%index D with a vector that describes the rows and clolumns we want 
 %isfinits is putting zeros in a new vector where there are non-finite values 
                        % and find is finding the values where there are ones and not zeros 

for x = 1:1:3
    for y= 1:1:3
    Iy1 = y:nx:nr;
    Ix1 = x:nx:nc;
    D_idx = D(Ix1, Iy1);
%     D_idx = D_idx(isfinite(D_idx));
%     D_idx = D_idx(:);
    Didx_mu_30(x,y) = nanmean(D_idx(:));
    Didx_min_30(x,y) = nanmin(D_idx(:));
    Didx_max_30(x,y) = nanmax(D_idx(:));
    Didx_std_30(x,y) = nanstd(D_idx(:));
    end 
end
D_mu = mean(Didx_mu_30, 'all');
D_std = std(Didx_mu_30(:));
x0 =(D_mu-100*D_std):D_std/100:(D_mu+100*D_std);
figure
% histogram(Didx_mu_30, 10,'Normalization','pdf'); hold on 
[S, nbins] = hist(Didx_std_30, 10)
dx = mean(diff(nbins))
N = length(Didx_std_30)
f = S*dx/N
bar(nbins,f)

gauss_c = mypdf(x0, D_mu, D_std);
plot(x0, gauss_c, 'r') 
%% Question 14
clear D,Ix, Ix1, Iy1, D_idx, Dstd, Dmax, Dmin, Dmean;
D = load('elevations.txt')

 %isfinits is putting zeros in a new vector where there are non-finite values 
                        % and find is finding the values where there are ones and not zeros 

nr = 111
nc = 110
samplesize = [50, 100, 200, 400,500]
for q= 1:length(samplesize)
    nx = round(samplesize(q)/9);

    for x = 1:1:nx
        for y = 1:1:nx
        Iy1 = y:nx:nr;
        Ix1 = x:nx:nc;
        D_idx = D(Ix1, Iy1);
%         D_idx = (isfinite(D_idx));
%         D_idx = D_idx(:);
        Didx_mu(x,y) = nanmean(D_idx(:));
        Didx_min(x,y) = nanmin(D_idx(:));
        Didx_max(x,y) = nanmax(D_idx(:));
        Didx_std(x,y) = nanstd(D_idx(:));
        end
     
    end 
    x_mu(q) = nanmean(Didx_mu(:));
    x_std(q) = nanstd(Didx_mu(:))
    len(q) = length(Didx_mu) % number of samples 
    
end
figure
plot(len, x_mu); hold on
plot(len, x_mu-x_std); hold on
plot(len, x_mu+x_std); 
% len, x0,x_std
%  Dmean = nanmean(Didx_mu, 'all'); % mean of D2
%  Dmin = nanmin(Didx_mu,[], 'all') ;
%  Dstd = nanstd(Didx_mu(:));
%  Dmax = nanmax(Didx_mu, [], 'all') ;

% figure(15)
% subplot(2,2,1)
% plot(len, x_mu)
% title('D Mean')
% subplot(2,2,2)
% boxplot(len, x_mu-x_std)
% title('D Standard Deviation')
% subplot(2,2,3)
% boxplot(len, x_mu+x_std)
% title('D Maximum')
% subplot(2,2,4)
% boxplot(Dmin)
% title('D Minimum')

%% GAUSSIAN FUNCTION 

function f =mypdf(z,mu,sig)

A = 1/(sig*sqrt(2*pi));
B = (z-mu).^2;
C = 2*sig.^2;
f = A*exp(-B./C);
end