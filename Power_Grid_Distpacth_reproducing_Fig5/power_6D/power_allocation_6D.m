
function cost=power_allocation_6D(x)% 
    penalty=10000;% 
    popMax=[500,200,300,150,200,120];
    popMin=[100,50,80,50,50,50];
    a=[0.007,0.0095,0.009,0.009,0.008,0.0075];
    b=[7,10,8.5,11,10.5,12];
    c=[240,200,220,200,220,190];
    e=[300,200,150,150,150,150];
    f=[0.031,0.042,0.063,0.063,0.063,0.063];
    cost=0;
    for i=1:6
        cost=cost+a(i)*x(i)*x(i)+b(i)*x(i)+c(i)+abs( e(i)*sin(f(i)*(popMin(i)-x(i))) );
    end
    if x(6)<popMin(6) || x(6)>popMax(6)
        cost=cost+penalty;
    end
end