
function object_value = kernel_gaussian(b, y, x, h,order)


% This function is to calculate the value of higher order Gaussian kernel
% h is the bandwidth
% order is the kernel order
% Reference: Hansen, B. E. (2005). Exact mean integrated squared error of higher order kernel estimators.

dum = (y - x*b)/h;

if order==6
    K=-1/8*(15-10*dum.^2+dum.^4).*exp(-dum.^2/2)/sqrt(2*pi);
elseif order==8
    K=-1/48*(105-105*dum.^2+21*dum.^4-dum.^6).*exp(-dum.^2/2)/sqrt(2*pi);
elseif order==10
    K=-1/384*(945-1260*dum.^2+378*dum.^4-36*dum.^6+dum.^8).*exp(-dum.^2/2)/sqrt(2*pi);
elseif order==12
    K=-1/3840*(10395-17325*dum.^2+6930*dum.^4-990*dum.^6+55*dum.^8-dum.^(10)).*exp(-dum.^2/2)/sqrt(2*pi);
end
object1=K/h;
object_value = sum(object1);

end