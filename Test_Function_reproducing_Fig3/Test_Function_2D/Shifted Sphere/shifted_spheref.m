function [y] = shifted_spheref(xx,oo)
xx = xx - oo;
d = length(xx);
sum = 0;
for ii = 1:d
	xi = xx(ii);
	sum = sum + xi^2;
end

y = sum;
end
