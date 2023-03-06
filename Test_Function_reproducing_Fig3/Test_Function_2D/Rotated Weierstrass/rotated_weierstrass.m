function [f]=rotated_weierstrass(x,M)
[ps,D]=size(x);
x = reshape(x,1,D);
x=x*M;
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c = -D*sum(c1 .* cos(c2.*0.5));
% sum(c1 .* cos(c2.*x(i)));
% c=-ww(0.5,c1,c2);
for i=1:D
    f=f+sum(c1 .* cos(c2.*x(i)));
end
f=f+c;


