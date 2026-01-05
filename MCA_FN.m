function [rank,ric]=MCA_FN(X,h,tol,rmax,rep)

% This function is to estimate the number of factor with modal component analysis
% (MCA) method.

% X is the T by X matrix of observed varibales
% h is the bandwidth
% tol is the tolerance for convergence
% rmax is the maximum number of factors
% rep is the number of different initializations

% rank is the estimated number factors using the rank estimator
% ric is the estimated number factors using the information-criteria estimator
% Reference: Sun, Z. and Tu, Y. (2025). MCA: High-Dimensional Modal Component Analysis towards the Mode

VALUE=[];
[T,N]=size(X);
L=min(T,N);

% estimate the factor and loadings with differnet factor number 
for q=1:rmax
    [Fhat,Lhat]=MM(X,q,h,tol,rep);
    value=mean(mean(kernel(X,Fhat*Lhat',h)));
    VALUE=[VALUE,value];
end

% calculate the IC
   [F0,L0]=MM(X,0,h,tol,1);
   value0=mean(mean(kernel(X,F0*L0',h)));
   array=1:rmax;
   IC=-VALUE+(VALUE(:,1)-value0)*L^(-4/7*1/2).*array;
   [temp,ric]=min(IC);
   Lrmax=diag(Lhat'*Lhat/N);
   thre=Lrmax(1)*(1/L)^(4/7*7/15); %threshold
   rank=sum(Lrmax>thre);

   
 
   
   