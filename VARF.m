function var_fhat=VARF(X, Fhat,Lhat,h,t)

% This function is to calculate the asymptotic variance of MCA factors
% X is the observation matrix
% Fhat,Lhat are the estimated MCA fcator and loading matrix
% h is the bandwidth
% t is the time index 
% Reference: Sun, Z. and Tu, Y. (2025). MCA: High-Dimensional Modal Component Analysis towards the Mode

[T,N]=size(X);

Ehat=(X(t,:)-Fhat(t,:)*Lhat')';  
E1=kernel2(Ehat,0,h);
E2=Lhat'*(E1.*Lhat)/N;
E3=kernel1(Ehat,0,h).^2/h;
E4=Lhat'*(E3.*Lhat)/N;
E5=pinv(E2)*E4*pinv(E2);    
var_fhat=E5;

end
