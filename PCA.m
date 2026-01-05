function [Fhat,Lhat]=PCA(X,r)

% This function is to estimate the model with principle component
% analysis (PCA) method
% X is the T by N matrix of observed variables 
% r is the number of factors
% Fhat is the T by r matrix of estimated PCA factors
% Lhat is the N by r matrix of estimated PCA factor loadings
%Reference: Bai, J. and Ng, S. (2002). Determining the number of factors in approximate factor models.

[T,N]=size(X);
[Y,W,V]=svd(X);
Fhat= Y(:,1:r)*sqrt(T);
Lhat= X' * Fhat/T;
[Fhat,Lhat]=normalize(Fhat,Lhat);