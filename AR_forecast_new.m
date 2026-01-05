function [y_pre,optimal_lag] = AR_forecast_new(Y,h,lag_max)

% This function is to generate h-step-ahead forecast of y based on autoregressive regression
% Y is a vector of univariate time series
% h is the forecast horizon
% lag_max is the maximum length of lags
% y_pre is the h-step-ahead forecast of y
% optimal_lag is the opimal lag order determined by BIC

T=length(Y); 

Y_h=Y(lag_max+h:T);

T_sample=T-h; % the size of data used to estimate the model
% 
% Y_h=[];
% for j=lag_max+1:T-h+1
%     Y_h=[Y_h; mean(Y(j:j+h-1))];
% end
% 
n=length(Y_h);

X=[];
for j=lag_max:-1:1
    X=[X Y(j:j+n-1)  ];
end
X=[ones(n,1) X];

% calculate the SSRs for each k=1:lag_max
SSR=[];
for k=1:lag_max
    dum_X = X(:,1:k+1);
    PX= dum_X* inv(dum_X'*dum_X)* (dum_X');
    Resid = (eye(n) - PX)*Y_h;
    SSR=[SSR; Resid'*Resid];
end

% use BIC to choose the optimal lag length
K=1:k;
BIC = log(SSR/n)+ (K'+1)*log(n)/n;

[dum, optimal_lag] = min(BIC);  


% generate the forecast
dumX=X(:,1:optimal_lag+1);
beta = inv(dumX'*dumX)* (dumX')*Y_h;

%T1=T-h+1;
T1=T;
y_pre = (beta')*[1;Y(T1:-1:T1-optimal_lag+1)];












