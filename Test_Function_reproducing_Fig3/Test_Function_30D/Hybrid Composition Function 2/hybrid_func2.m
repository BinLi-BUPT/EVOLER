function fit=hybrid_func2(x)
persistent  fun_num func o sigma lamda
    x = reshape(x,1,30); 
    [ps,D]=size(x);
    fun_num = 10; 
    o = zeros(1,D);
    sigma=ones(1,fun_num);
    bias=((1:fun_num)-1).*0;
    func.f1=str2func('fsphere');  
    func.f2=str2func('fsphere'); 
    func.f3=str2func('fweierstrass'); 
    func.f4=str2func('fweierstrass');
    func.f5=str2func('frastrigin');
    func.f6=str2func('frastrigin');
    func.f7=str2func('fgriewank');
    func.f8=str2func('fgriewank');
    func.f9=str2func('flevy');
    func.f10=str2func('flevy');
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
%     for i=1:fun_num
%         eval(['M.M' int2str(i) '=diag(ones(1,D));']);
%     end
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda);
end


function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda)
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
tmp1 = 100*fsphere((x-o)./repmat(lamda(1,:),ps,1))/fsphere(x1./lamda(1,:));
tmp2 = 100*fweierstrass((x-o)./repmat(lamda(3,:),ps,1))/fweierstrass(x1./lamda(3,:));
tmp3 = 100*frastrigin((x-o)./repmat(lamda(5,:),ps,1))/frastrigin(x1./lamda(5,:));
tmp4 = 100*fgriewank((x-o)./repmat(lamda(7,:),ps,1))/fgriewank(x1./lamda(7,:));
tmp5 = 100*flevy((x-o)./repmat(lamda(9,:),ps,1))/flevy(x1./lamda(9,:));
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
x = x + repmat(1, 1, D);
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

function f=flevy(xx)
d = length(xx);
xx = xx + repmat(1, 1, d);
for ii = 1:d
	w(ii) = 1 + (xx(ii) - 1)/4;
end

term1 = (sin(pi*w(1)))^2;
term3 = (w(d)-1)^2 * (1+(sin(2*pi*w(d)))^2);

sum = 0;
for ii = 1:(d-1)
	wi = w(ii);
        new = (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
	sum = sum + new;
end
f = term1 + sum + term3;
end
