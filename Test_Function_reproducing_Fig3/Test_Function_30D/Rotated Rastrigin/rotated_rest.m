function y=rotated_rest(x,M)
[~,D]=size(x);
x = reshape(x,1,D);
x=x*M;
d = length(x);
sum = 0;
for ii = 1:d
	xi = x(ii);
	sum = sum + (xi^2 - 10*cos(2*pi*xi));
end

y = 10*d + sum;
end
    
  