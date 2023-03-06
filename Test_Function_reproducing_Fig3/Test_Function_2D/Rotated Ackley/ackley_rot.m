% ackley rotated function
function f=ackley_rot(x,M)
[ps,D]=size(x);
x = reshape(x,1,D);
c=100;
x=x*M;
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);