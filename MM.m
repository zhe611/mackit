function [Fhat,Lhat] = MM(X,r,h,tol,rep)

% This function is to estimate factor model with modal component analysis
% (MCA) method, using the minorization-maximization (MM) algorithm and
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

bu=1;
cu=1;

% choose the starting value for the AMEM algorithm 
rng(2024); 

obj=[];
F=[];
L=[];


for s=1:rep 

rng(2024+s);     

rng("shuffle")

F1=rand(T,r);
L1=rand(N,r);

obj1=-1;

obj0=0;

g=0;

while g<=100 & (obj0-obj1)>tol

obj1=obj0;

CC1=F1*L1';

temp=-(X-CC1).^2/2/h^2;

weight=exp(temp);

X1=CC1-weight.*(-(X-CC1)/h)/bu/cu*h;

[F1,L1]=PCA(X1,r);

temp=-(X-F1*L1').^2/2/h^2;

obj0=sum(sum(exp(temp)))/N/T;

obj0-obj1;

g=g+1;

end

obj=[obj,obj1];
F=[F,F1];
L=[L,L1];

end

[value,p]=max(obj);

F0=F(:,(r*p-r+1):r*p);
L0=L(:,(r*p-r+1):r*p);

[Fhat,Lhat]=normalize(F0,L0);


end

    
 
    