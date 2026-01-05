
function rank=QFA_FN(X,rmax,tol,tau,rep)

% This is function is to calculate the number of QFA factors using the rank-minimization estimator

% X is the T by N matrix of observed varibales
% rmax is the maximum number of factors
% tol is the tolerance of convergence
% tau is the quantile level
% rank is the estimated number of quantile factors at the tau-th quantile

% Reference: Chen, L., Dolado, J. J., and Gonzalo, J. (2021). Quantile factor models.

[T,N]=size(X);
L=min(T,N);

% estimating the qunatile factors and loadings
 [Fhat,Lhat]=IQR(X,rmax,tol,tau,rep);
 Lrmax=diag(Lhat'*Lhat/N);
 thre=Lrmax(1)*L^(-1/3); %threshold
 rank=sum(Lrmax>thre);
   