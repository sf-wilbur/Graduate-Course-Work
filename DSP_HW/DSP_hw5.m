% lets do the synthesis equation! (equation 8-2)

x = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0];
figure(1); clf

% X = fft(x);
N = length(x);
K = length(x)/2+1;
ks = 0:K-1; % frequency vector... goes up to Nyquist
 
subplot(2,2,[1 2])
stem(0:N-1,x) % time series (non-repeating)
title('time series')
 
subplot(2,2,3)
plot(0:N-1,abs(X))
title('amplitude')
subplot(2,2,4)
plot(0:N-1,angle(X),'o')
title('phase')
ylim([-pi pi])

function X = dft(x)
% compute the DFT of a data vector x using basis functions
if size(x,1)<size(x,2)
x = x'; % transpose matrix if needed
end
K = length(x)/2 + 1; % number of frequency points to compute
fax = linspace(0,pi,K); % frequency axis in radians
n = (0:length(x)-1)'; % vector of time series indices (starting with 0)
for k=1:K % loop through each frequency (in radians)
ff = fax(k); % incremental frequency
c = cos(ff*n); % cosine basis function
s = sin(ff*n); % sine basis function
re(k) = sum(x.*c(:,k)); % real coefficient
im(k) = sum(x.*s(:,k)); % imaginary coefficient
end
X = re - im*(i); % complex spectrum
end 


