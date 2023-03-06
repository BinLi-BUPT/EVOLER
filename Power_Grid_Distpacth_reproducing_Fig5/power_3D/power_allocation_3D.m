function cost=power_allocation_3D(x)
    penalty=10000;
    popMax=[600,200,400];
    popMin=[100,50,100];
    a=[0.001562,0.004820,0.001940];
    b=[7.92,7.97,7.85];
    c=[561,78,310];
    e=[300,150,200];
    f=[0.0315,0.063,0.042];
    cost=0;
    for i=1:3
        cost=cost+a(i)*x(i)*x(i)+b(i)*x(i)+c(i)+abs( e(i)*sin(f(i)*(popMin(i)-x(i))) );
    end
    if x(3)<=popMin(3) || x(3)>=popMax(3)
        cost=cost+penalty;
    end
end