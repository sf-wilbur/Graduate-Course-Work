D = load ('icevelocity.txt')
T = D(:,1) 
V = D(:,2)
A= [min(V):1:max(V)]

for N = 1:length(A);
       Vm= A(N);
       R_err(N) = sqrt(mean((V-Vm).^2))
end 
plot(A,R_err)