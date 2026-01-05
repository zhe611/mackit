
function x=kernel(t1,t2,h)

% calculate the value of the kernel function evaluated at the point t1-t2

x=exp((t1-t2).^2./(-2*h^2))./(sqrt(2*pi)*h);

end






