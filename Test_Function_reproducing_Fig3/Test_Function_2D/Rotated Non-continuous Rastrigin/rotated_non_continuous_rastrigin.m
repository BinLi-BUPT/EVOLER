%Noncontinuous Rastraigin
function val= rotated_non_continuous_rastrigin(x,M)
[~,D]=size(x);
x = reshape(x,1,D);
x=x*M;
val=0;
for i=1:D
    if abs(x(i))<1/2
        y=x(i);
    else
        y=round(2*x(i))/2;
    end
    val=val+y^2-10*cos(2*pi*y)+10;
end
end