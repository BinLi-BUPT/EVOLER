function [f]=weierstrass(x)
[ps,D]=size(x);
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c = -D*sum(c1 .* cos(c2.*0.5));
for i=1:D
    f=f+sum(c1 .* cos(c2.*x(i)));
end
f=f+c;


