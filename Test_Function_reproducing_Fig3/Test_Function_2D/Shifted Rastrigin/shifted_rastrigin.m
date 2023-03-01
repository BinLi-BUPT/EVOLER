function [y] = shifted_rastrigin(xx,oo)
xx = xx - oo;
d = length(xx);
sum = 0;
for ii = 1:d
	xi = xx(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end

y = 10*d + sum;

end

