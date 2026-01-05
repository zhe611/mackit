function [Fhat,Lhat] = IQR(X,r,tolerate,tau,rep)

% This function is to estimate the model using the quantile factor analyis (QFA) method.

% X is the T by N matrix of observed variables 
% r is the number of factors
% tolerate is the convergence threshold 
% tau is the quantile level
% rep is the number of intilizations
% Fhat is the T by r matrix of estimated QFA factors
% Lhat is the N by r matrix of estimated QFA loadings

% Reference: Chen, L., Dolado, J. J., and Gonzalo, J. (2021). Quantile factor models.

[T,N]=size(X);
obj=[];
F=[];
L=[];



for s=1:rep 
rng(2024+s);  

obj_1=1; 
obj_0=0;
step=0;

F0=rand(T,r)'; % intial value of the factor

while abs(obj_1-obj_0)>tolerate & step<=100
    lambda0=[];
    for i=1:N
        lambda0=[lambda0 rq_fnm(F0',X(:,i),tau)];
    end;
    obj_0 = mean( mean ( ( tau - (X< F0'*lambda0 )).*(X-F0'*lambda0))) ; 
        
    F1=[];
    for j=1:T
        F1=[F1 rq_fnm(lambda0',X(j,:)',tau)];
    end;
    F0=F1;
    
    lambda1=[];
    for i=1:N
        lambda1=[lambda1 rq_fnm(F1',X(:,i),tau)];
    end;
    obj_1 = mean( mean ( ( tau - (X< F1'*lambda1 )).*(X-F1'*lambda1))) ; 
   
step=step+1;
end

Fhat=F1';
Lambdahat=lambda1';

obj=[obj,obj_1];
F=[F,Fhat];
L=[L,Lambdahat];



end

[value,p]=min(obj);

Fhat=F(:,(r*p-r+1):r*p);
Lambdahat=L(:,(r*p-r+1):r*p);

% normalize the factors and loadings 
[Fhat,Lhat]=normalize(Fhat,Lambdahat);

