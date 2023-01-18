function [h,V,c, Iu, Iv] = semivariogram(x,y)
% simple vairogram function for equally spaced data 
% Input 
% x = distance vector 
% y = measurement vecotr 
% Output; 
% h = lag distance 
% y=semivariagram result 
% SNTX: [h,y] = semivariogram(x,y)
dx = mean(diff(x)); %average spacing 
extent = (max(x) - min(x)); % extent 
N = length(x); %number of data points 
h = dx:dx:extent/2 ; % lags --> only calculating to half the extent to avoid bias 
npairs = zeros(length(h), 1); %preallocate npairs 
V = zeros(length(h),1); %preallocate semivariance 
muz = mean(y); 
for n = 1:length(h) %loop over lags 
    npairs(n) = N-n; % number of pairs at each lag 
    Iu = 1:(N-n); % define the index of the head
    Iv = (n+1):N; % index to tails 
    V(n) = 1./(2*npairs(n))*sum((y(Iu)-y(Iv)).^2); %semivariance 
    c(n) = 1./(npairs(n)-1)*sum((y(Iu)-muz).*(y(Iv)-muz)); %covariance 

end 