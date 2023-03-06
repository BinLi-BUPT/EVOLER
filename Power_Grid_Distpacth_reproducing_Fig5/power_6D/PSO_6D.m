clear all;
clc;
format long;
runs=50;
c1=1.4962;             
c2=1.4962;             
maxGen=500;            
D=5;                  
N=20;                  
popSum=1260;
popMax=[500,200,300,150,200];
popMin=[100,50,80,50,50];
DMin=50; DMax=120;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;
for m=1:runs
    m
    for i=1:N
        while(1)
            pop(i,:)=popMin+(popMax-popMin).*rand(1,D);  
            x=pop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            if x(D+1)<=DMax && x(D+1)>=DMin, break; end
        end
        fitness(i)=power_allocation_6D(x);
        V(i,:)=vMin+(vMax-vMin).*rand(1,D); 
    end
    %----------------------------
    [fitnessgbest bestindex]=min(fitness);
    yy(1)=fitnessgbest;
    gbest=pop(bestindex,:);
    pbest=pop;
    fitnesspbest=fitness;
    w=0.7968; 
    for i=1:maxGen
        for j=1:N
            V(j,:)=w*V(j,:)+c1*rand(1,D).*(pbest(j,:)-pop(j,:))+c2*rand(1,D).*(gbest-pop(j,:));
            for d=1:D
                if V(j,d)>vMax(d), V(j,d)=vMax(d);end
                if V(j,d)<vMin(d), V(j,d)=vMin(d);end
            end
            pop(j,:)=pop(j,:)+V(j,:);
            for d=1:D
                if pop(j,d)>popMax(d), pop(j,d)=popMax(d); end
                if pop(j,d)<popMin(d), pop(j,d)=popMin(d); end
            end
            x=pop(j,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            fitness(j)=power_allocation_6D(x);


            if fitness(j)<fitnesspbest(j)
                pbest(j,:)=pop(j,:);
                fitnesspbest(j)=fitness(j);
            end

            if fitness(j)<fitnessgbest
                gbest=pop(j,:);
                fitnessgbest=fitness(j);
            end

        end
        yy(i+1)=fitnessgbest;  
    end
%     plot(yy,'r');
%     hold on;
%     xlabel('Generation');
%     ylabel('fitness');
%     title('Power Allocation PSO-6D ');
    COST(m,1)=fitnessgbest;
    COST(m,1)
    COSTS(m,:)= yy;
end
save('PSO_6D','COSTS')


