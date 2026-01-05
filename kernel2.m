function x=kernel2(t1,t2,h)

% calculate the second deriative of the kernel function evaluated at the point t1-t2

t=t1-t2;
x=h^(-3)*(t.^2./h^2-1)./sqrt(2*pi).*exp(-t.^2/2/h^2);
end