function x=kernel1(t1,t2,h)

%calculate the first derivative of the kernel function

t=t1-t2;
x=-t./h./sqrt(2*pi).*exp(-t.^2/2/h^2);

end