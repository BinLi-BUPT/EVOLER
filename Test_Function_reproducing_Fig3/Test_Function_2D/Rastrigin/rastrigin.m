function y = rastrigin(x)
%     x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
%     f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end

y = 10*d + sum;

end