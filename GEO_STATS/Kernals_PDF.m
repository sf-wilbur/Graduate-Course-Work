%hands on kernal PDF
D = load('elevations.txt')
D = D(isfinite(D)); %get rid of NANs
[nc,xbins] = hist(D,100); % calculate histogram 
dx = xbins(2)-xbins(1)%bin width 
frdh = nc/(sum(nc*dx)); %relative density 
figure(1); clf
bar(xbins,frdh);

dx = 0.5 %resolution of the pdf
xm = min(xbins):dx:max(xbins);


h = 10; 
for n=1:length(xm)
    x0 = xm(n);
    dist = D-x0; %distance from x0 to all other data 
    Ix = find(abs(dist)<h) % finding data within window h
    w =15/16*(1-dist(Ix)/h).^2
    kde(n) = sum(w); %store the estiamte for each position in our array 
    
end 

kde= kde/(sum(kde*dx)) % this is the normalization of kde 
figure(1); hold on 
plot(xm, kde, 'r', 'linewidth', 2)

kernal_pdf(D,xm,h)

