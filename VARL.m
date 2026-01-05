function var_Lhat=VARL(X, Fhat,Lhat,h,i)

% This function is to calculate the asymptotic variance of MCA loadings
% X is the observation matrix
% Fhat,Lhat are the estimated MCA fcator and loading matrix
% h is the bandwidth
% i is the unit index 
% Reference: Sun, Z. and Tu, Y. (2025). MCA: High-Dimensional Modal Component Analysis towards the Mode

[T,N]=size(X);

Ehat=(X(:,i)-Fhat*Lhat(i,:)');  
E1=kernel2(Ehat,0,h);
E2=Fhat'*(E1.*Fhat)/T;
E3=kernel1(Ehat,0,h).^2/h;
E4=Fhat'*(E3.*Fhat)/T;
E5=pinv(E2)*E4*pinv(E2);    
var_Lhat=E5;

end
