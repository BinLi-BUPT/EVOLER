function fit=shifted_hybrid_func1(x,oo)
persistent  fun_num func o sigma lamda
x = reshape(x,1,30);
[ps,D]=size(x);
fun_num = 10;
o = oo;
sigma=ones(1,fun_num);
lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
lamda=repmat(lamda,1,D);

fit=hybrid_composition_func(x,fun_num,o,sigma,lamda);
end


function fit=hybrid_composition_func(x,fun_num,o,sigma,lamda)
[ps,D]=size(x);
for i=1:fun_num
    %     oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-o).^2,2)./2./(D*sigma(i)^2));
end

% [tmp,tmpid]=sort(weight,2);
% for i=1:ps
%     weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
% end
weight=weight./repmat(sum(weight,2),1,fun_num);
fit=0;
x1=5*ones(1,D);
tmp1 = 2000*frastrigin((x-o)./repmat(lamda(1,:),ps,1))/frastrigin(x1./lamda(1,:));
tmp2 = 2000*fweierstrass((x-o)./repmat(lamda(3,:),ps,1))/fweierstrass(x1./lamda(3,:));
tmp3 = 2000*fgriewank((x-o)./repmat(lamda(5,:),ps,1))/fgriewank(x1./lamda(5,:));
tmp4 = 2000*felliptic((x-o)./repmat(lamda(7,:),ps,1))/felliptic(x1./lamda(7,:));
tmp5 = 2000*fsphere((x-o)./repmat(lamda(9,:),ps,1))/fsphere(x1./lamda(9,:));
fit = (weight(:,1)*tmp1+weight(:,2)*tmp1+weight(:,3)*tmp2+weight(:,4)*tmp2+...
    weight(:,5)*tmp3+weight(:,6)*tmp3+weight(:,7)*tmp4+weight(:,8)*tmp4+...
    weight(:,9)*tmp5+weight(:,10)*tmp5);


% for i=1:fun_num
%     oo=repmat(o(i,:),ps,1);
%     eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1));']);
%     eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))']);
%     fit1=2000.*f./f1;
%     fit=fit+weight(:,i).*(fit1+bias(i));
% end
end
function f=fsphere(x)
[ps,D]=size(x);
f=sum(x.^2,2);
end
%--------------------------------
function f=fgriewank(xx)
[ps,D]=size(xx);
f=1;
for i=1:D
    f=f.*cos(xx(:,i)./sqrt(i));
end
f=sum(xx.^2,2)./4000-f+1;
end
%--------------------------------
function f=frastrigin(x)
d = length(x);
sum = 0;
for ii = 1:d
    xi = x(ii);
    sum = sum + (xi^2 - 10*cos(2*pi*xi));
end
f = 10*d + sum;
end
%--------------------------------
function [f]=fweierstrass(x)
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
end

function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
    f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
end
