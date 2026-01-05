function [Fhat,Lhat] = AM_Gaussian(X,r,h,tol,rep)

% This function is to estimate factor model with modal component analysis
% (MCA) method, using the alternating maximization (AM) algorithm and
% gaussian kernel.

% X is the T by N matrix of observed variables 
% r is the number of factors
% h is the bandwidth 
% tol is the convergence threshold 
% rep is the number of different initializations
% Fhat is the T by r matrix of estimated MFA factors
% Lhat is the N by r matrix of estimated MFA factor loadings
% Reference: Sun, Z. and Tu, Y. (2025). MCA: High-Dimensional Modal Component Analysis towards the Mode


[T,N]=size(X);

% choose the starting value for the AMEM algorithm 
rng(2024); 
obj=[];
F=[];
L=[];

for s=1:rep 

rng(2024+s);     

F0=rand(T,r);
Lamb0=rand(N,r);

obj_1=1; 
obj_0=0;
t=[];
t1=[];

step=0; %caluate the iteration time

while abs(obj_1-obj_0)>tol & step<=100
    for i=1:N
        Lamb0(i,:)=(mdreg(F0,X(:,i),Lamb0(i,:)',h,tol))';
    end;
    obj_0=mean(mean(kernel(X-F0 * Lamb0',0,h)));

    F1=normrnd(0,1,[T,r]);
    for j=1:T
        F1(j,:)=(mdreg(Lamb0,(X(j,:))',(F0(j,:))',h,tol))';
    end;
    
    lamb1=normrnd(0,1,[N,r]);
    for i=1:N
        Lamb1(i,:)=(mdreg(F1,X(:,i),Lamb0(i,:)',h,tol))';
    end;
    obj_1 = mean(mean(kernel(X-F1 * Lamb1',0,h)));
    

    F0=F1;
    Lamb0=Lamb1;
    step=step+1;
end
obj=[obj,obj_1];
F=[F,F1];
L=[L,Lamb1];
end

% find the best local maxima
[value,p]=max(obj);

F0=F(:,(r*p-r+1):r*p);
L0=L(:,(r*p-r+1):r*p);


% normalize the factors and loadings
[Fhat,Lhat]=normalize(F0,L0);

return;

end

    
 
    