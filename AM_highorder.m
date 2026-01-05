function [Fhat,Lhat,step] = AM_highorder(X,r,h,order,tol,rep)

% This function is to estimate factor model with modal component analysis
% (MCA) method, using the alternating maximization (AM) algorithm and
% higher-order gaussian kernel.

% X is the T by N matrix of observed variables 
% r is the number of factors
% h is the bandwidth 
% order is the order of the kernel function
% tol is the convergence threshold 
% rep is the number of different initializations
% Fhat is the T by r matrix of estimated MFA factors
% Lhat is the N by r matrix of estimated MFA factor loadings
% Reference: Sun, Z. and Tu, Y. (2025). MCA: High-Dimensional Modal Component Analysis towards the Mode

options = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton');

[T,N]=size(X);

% choose the starting value for the AMEM algorithm 
rng(2024); 
obj=[];
F=[];
L=[];
step=0;

for s=1

rng(2024+s);     

%F0=rand(T,r);
%L0=rand(N,r);

[F0,L0]=MM(X,r,h,tol,rep); % use the MM estimates as the initial value

obj_1=1; 
obj_0=0;
t=[];
t1=[];
step=0; % calculate the iteration time

while abs(obj_1-obj_0)>tol & step<=100
    
    L1=[];
    
    for i=1:N
    obj_fun = @(b) kernel_gaussian(b, X(:,i), F0, h, order);
    b0 = L0(i,:)';
    [b_opt, fval] = fminunc(obj_fun, b0,options);
    L1=[L1;b_opt'];
    end

    obj_0=mean(mean(kernel(X-F0 * L1',0,h)));


    F1=[];
    for j=1:T
    obj_fun = @(b) kernel_gaussian(b, X(j,:)', L1, h, order);
    b0 = F0(j,:)';
    [b_opt, fval] = fminunc(obj_fun, b0,options);
    F1=[F1;b_opt'];
    end

    obj_1 = mean(mean(kernel(X-F1 * L1',0,h)));
    
    F0=F1;
    L0=L1;

    step=step+1;
end
obj=[obj,obj_1];
F=[F,F1];
L=[L,L1];

end

% find the best local maxima
[value,p]=max(obj);

F0=F(:,(r*p-r+1):r*p);
L0=L(:,(r*p-r+1):r*p);

% normalize the factors and loadings
[Fhat,Lhat]=normalize(F0,L0);

return;

end

    
 
    