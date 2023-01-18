function [h,V,npairs] = semivariogram_mc(x,y,np)
% simple vairogram function for equally spaced data 
% Input 
% x = distance vector 
% y = measurement vecotr 
% Output; 
% h = lag distance 
% y=semivariagram result 
% SNTX: [h,y, npairs] = semivariogram(x,y, np)
dx = mean(diff(x)); %average spacing 
extent = (max(x) - min(x)); % extent 
N = length(x); %number of data points 
h = dx:dx:extent/2 ; % lags --> only calculating to half the extent to avoid bias 
npairs = zeros(length(h), 1); %preallocate npairs 
V = zeros(length(h),3); %preallocate semivariance 
for n = 1:length(h) %loop over lags 
    npairs(n) = N-n; % number of pairs at each lag 
    Iu = 1:(N-n); % define the index of the head
    Iv = (n+1):N; % index to tails 
    V(n) = 1/(2*npairs(n))*sum((y(Iu)-y(Iv)).^2); %semivariance 
    Vt = zeros(10,1); 
    for m =1:10 %monte carlo for uncertainties 
        I2 = randsample(Iu, np); %randsample of pairs 
        Iut = Iu(I2); 
        Ivt = Iv(I2); 
        Vt(m) = 1/(2*npairs(n))*sum((y(Iut)-y(Ivt)).^2); % semivariance of randsample
    end 
    V(n,:) = quantile(Vt,[0.025 0.5 0.975]); 
end 