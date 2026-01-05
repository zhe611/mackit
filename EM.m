
function  Z= EM(X,y,beta,h)

% This function is to conduct E-M iterations in modal regression
% X is the matrix of regressors
% y is a vector of observed variables
% beta a vector of regression coefficient
% h is the bandwidth
% Z is a vector of updated regression coefficient
% Reference: Yao, W. and Li, L. (2014). A new regression model: modal linear regression.


% calculate the weight
weight=kernel(y,X*beta,h);
W = diag( weight );
temp= X' * W * X;
P=pinv(temp);

% update the value
Z= P * X' * W * y;

end








