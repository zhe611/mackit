
function ic=PCA_FN(X,rmax)

% This function is to estimate the number PCA factors using the infromation criteria mrthod
% X is a T by N observed variables
% rmax is the maximum number of factors
%Reference: Bai, J. and Ng, S. (2002). Determining the number of factors in approximate factor models.


[T,N]=size(X);
NT=N*T;
NT1=N+T;
L=min(N,T);
PC=zeros(1,rmax); %kmax+1列
Sigma=zeros(1,rmax);

for i=1:rmax
[Fhat,lambda]=PCA(X,i);    
ehat=X-Fhat*lambda'; %error
Sigma(i)=mean(mean(ehat.^2));
end
array=1:rmax;

PC=Sigma+log(NT/NT1)*NT1/NT*Sigma(rmax).*array; %penalty


[temp,ic]=min(PC);


