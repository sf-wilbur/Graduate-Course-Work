clear all; close all; 

%% Lecture 11/04/2021 Class Example 
load('SMPstability.mat');
WLS = X(:,34); %strength of the layer where failure occured 
I2 = find(S ==1); %index to unstable conditions 
I3 = find(S==0); %stable indicies 
figure(1); clf 
boxplot(WLS, S, 'notch', 'on') 

xmod=min(WLS):0.001:max(WLS);  %places where we are going to evaluate 
funstable= ksdensity(WLS(I2),xmod, 'width', 0.01); %conditional probabilty 
fstable= ksdensity(WLS(I3),xmod, 'width', 0.01); %stable probality 
Pus = length(I2)/length(WLS); %percent of unstable 
Ps = length(I3)/length(WLS); %percent of stable 
figure(2); clf 
plot(xmod, Pus*funstable);hold on
plot(xmod, Ps*fstable, 'r')
% test different cut values 
xc = 0:0.1:0.5;
Sm= zeros(size(S)); %set model estimates to be the same size of my observed values for comparison 
for n = 1:length(xc) 
    I4 =find(WLS<xc(n)); %find values we will classify as unstable 
    I5 =find(WLS>=xc(n)); % find values we will classify as stable 
    Sm(I4) = 1; %set modeled values as unstable 
    TP = find(Sm==1 & S==1);
    FP = find(Sm ==1 & S ==0);
    TN = find(Sm==0 & S==0);
    FN = find(Sm==0 & S ==1);
    TA(n) = (length(TP)+length(TN))./length(S); %total accuracy 
    SA(n) = length(TN)/length(I3);
    UA(n) = length(TP)/length(I2); %unstable accuracy 

end 
figure(4); clf 
plot(xc,TA); hold on; plot(xc,SA, 'r'); plot(xc, UA, 'g')
legend('TA', 'SA', 'US');
% new code that uses fitree.m instead of classregtree.m
t2 = fitctree(X,S,'prior', 'empirical', 'PredictorNames', Xnames)
view(t2, 'mode', 'graph');
m = max(t2.PruneList)-1; %testing all levels except full tree 
[E, ~,~,bestLevel]=cvloss(t2, 'SubTrees',0:m, 'KFold',10);
t2prune = prune(t2, 'Level', bestLevel); %prune the tree to the optimal level based on CV 
view(t2prune, 'Mode', 'graph'); 
