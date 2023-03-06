function f=rotated_griewank(x,M)
[~,D]=size(x);
x = reshape(x,1,D);
x=x*M;
f1=1;
for i=1:D
    f1=f1.*cos(x(:,i)./sqrt(i));
end
f2=0;
for i=1:D
    f2=f2+x(:,i)^2;
end
f2=f2/4000;
f=f2-f1+1;
end

