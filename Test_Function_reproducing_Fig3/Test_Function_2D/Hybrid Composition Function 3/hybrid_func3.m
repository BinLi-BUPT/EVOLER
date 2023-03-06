function fit=hybrid_func3(x)
persistent  fun_num func o sigma lamda bias
x = reshape(x,1,2);
[ps,D]=size(x);
fun_num = 10;
load hybrid_func1_data % saved the predefined optima
if length(o(1,:))>=D
    o=o(:,1:D);
else
    o=-5+10*rand(fun_num,D);
end
bias=((1:fun_num)-1).*100;
sigma=ones(1,fun_num);
lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
lamda=repmat(lamda,1,D);
fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias);
end


function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias)
[ps,D]=size(x);
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2));
end

[tmp,tmpid]=sort(weight,2);
for i=1:ps
    weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
end
weight=weight./repmat(sum(weight,2),1,fun_num);
fit=0;
x1=5*ones(1,D);
for i=1:fun_num
    oo(i,:)=repmat(o(i,:),ps,1);
end

tmp1 = 2000*frastrigin((x-oo(1,:))./repmat(lamda(1,:),ps,1))/frastrigin(x1./lamda(1,:));
tmp2 = 2000*frastrigin((x-oo(2,:))./repmat(lamda(2,:),ps,1))/frastrigin(x1./lamda(2,:));

tmp3 = 2000*fweierstrass((x-oo(3,:))./repmat(lamda(3,:),ps,1))/fweierstrass(x1./lamda(3,:));
tmp4 = 2000*fweierstrass((x-oo(4,:))./repmat(lamda(4,:),ps,1))/fweierstrass(x1./lamda(4,:));
tmp5 = 2000*fgriewank((x-oo(5,:))./repmat(lamda(5,:),ps,1))/fgriewank(x1./lamda(5,:));
tmp6 = 2000*fgriewank((x-oo(6,:))./repmat(lamda(6,:),ps,1))/fgriewank(x1./lamda(6,:));
tmp7 = 2000*fackley((x-oo(7,:))./repmat(lamda(7,:),ps,1))/fackley(x1./lamda(7,:));
tmp8 = 2000*fackley((x-oo(8,:))./repmat(lamda(8,:),ps,1))/fackley(x1./lamda(8,:));
tmp9 = 2000*fsphere((x-oo(9,:))./repmat(lamda(9,:),ps,1))/fsphere(x1./lamda(9,:));
tmp10 = 2000*fsphere((x-oo(10,:))./repmat(lamda(10,:),ps,1))/fsphere(x1./lamda(10,:));
fit = (weight(:,1)*(tmp1+bias(1))+weight(:,2)*(tmp2+bias(2))+weight(:,3)*(tmp3+bias(3))+weight(:,4)*(tmp4+bias(4))+...
    weight(:,5)*(tmp5+bias(5))+weight(:,6)*(tmp6+bias(6))+weight(:,7)*(tmp7+bias(7))+weight(:,8)*(tmp8+bias(8))+...
    weight(:,9)*(tmp9+bias(9))+weight(:,10)*(tmp10+bias(10)));

end

function f=fsphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D]=size(x);
f=sum(x.^2,2);
end
%--------------------------------
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
%--------------------------------
end
function f=fgriewank(x)
[ps,D]=size(x);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;

end
%--------------------------------
function f=fackley(x)
[ps,D]=size(x);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end
%--------------------------------
function f=frastrigin(x)
[ps,D]=size(x);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end
%--------------------------------
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
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
c=-w(0.5,c1,c2);
for i=1:D
    f=f+w(x(:,i)',c1,c2);
end
f=f+c*D;
end
function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
    y(k) = sum(c1 .* cos(c2.*x(:,k)));
end
end
%--------------------------------
function f=fE_ScafferF6(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
f=f+feval(fhd,x(:,[D,1]));
end
%--------------------------------
function f=fE_ScafferF6_noncont(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
f=f+feval(fhd,x(:,[D,1]));
end
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
f=f+F8F2(x(:,[D,1]));
end

%--------------------------------
function f=fschwefel_102(x)
[ps,D]=size(x);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
end
%--------------------------------
function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
    f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
end
%--------------------------------
% classical Gram Schmid
function [q,r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid
%
[n,m] = size(A);
q = A;
for j=1:m
    for i=1:j-1
        r(i,j) = q(:,j)'*q(:,i);
    end
    for i=1:j-1
        q(:,j) = q(:,j) -  r(i,j)*q(:,i);
    end
    t =  norm(q(:,j),2 ) ;
    q(:,j) = q(:,j) / t ;
    r(j,j) = t  ;
end
end

function M=rot_matrix(D,c)
A=normrnd(0,1,D,D);
P=cGram_Schmidt(A);
A=normrnd(0,1,D,D);
Q=cGram_Schmidt(A);
u=rand(1,D);
D=c.^((u-min(u))./(max(u)-min(u)));
D=diag(D);
M=P*D*Q;
end
