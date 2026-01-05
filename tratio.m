
function x=tratio(F,F0)

% This function is to calculate the trace-ratio statistics of the matrix estimation
% F0 is the T by r matrix of ture matrix
% F is the T by r matrix of estimated matrix
% x is the value of the trace-ratio statistics of F

a=sum(diag( F0' * F0));
b=sum(diag( F0'*F*pinv(F'*F )*F'*F0));
x=b/a;

end





