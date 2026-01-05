
function beta2=mdreg(X,y,beta0,h,tol)

% This function is to conduct modal regression. 
% X is the matrix of regressors
% y is a vector of observed variables
% beta0 is a vector of initial value of the regression coeffcient
% h is the bandwidth
% tol is the convergence threshold 
% beta2 is a vector of updated regression coeffcient

% Reference: Yao, W. and Li, L. (2014). A new regression model: modal linear regression.


beta1=EM(X,y,beta0,h);
gap=beta1-beta0;
e=sqrt(gap'* gap);

% iterate until convergence
while e>tol
    beta0=beta1;
    beta1=EM(X,y,beta0,h);
    gap=beta1-beta0;
    e=sqrt(gap'* gap);
end
beta2=beta1;




