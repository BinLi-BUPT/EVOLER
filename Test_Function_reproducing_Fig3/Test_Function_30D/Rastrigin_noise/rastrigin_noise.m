function f = rastrigin_noise(x)
% [ps,D]=size(x);
% f=sum(x.^2,2).*(1+ (normrnd(0,1e-6,ps,1)));
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end
f = (10*d + sum)*(1+ 1e-7*(randn));