function rmse=model_variogram(h,V,c,a,type)
% variogram using the bounded linear and spherical models
% INPUT: h = lags at which modeled estimates are made
%        c = variogram sill = sig^2 (model parameter)
%        a = variogram range (model parameter)
%        V= experimental variogram estimate
%      type = 'L' for linear, 'S' for spherical

% first calculate variogram for lags less than range
ind=find(h<=a);
switch type 
    case'L'
        Vm(ind)=c*h(ind)/a; % bounded linear
    case 'S'
        Vm(ind)=c*(3*h(ind)/(2*a)-1/2*(h(ind)/a).^3); % spherical
end

% now define variogram for lags greater than range
ind2=h>a; % find points greater than range
Vm(ind2)=c; % set equal to sill
% V=V(:); 
% Vm=Vm(:); 
rmse = sqrt(mean((Vm(:)-V(:)).^2)); %root mean squared error
