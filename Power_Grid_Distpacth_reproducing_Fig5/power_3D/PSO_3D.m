clear all;
clc;
format long;
runs=50;  
maxGen=1000; 
c1=1.4962;
c2=1.4962;          
D=2;             
N=20;                 
popSum=850;
popMax=[600,200];
popMin=[100,50];
DMin=100;
DMax=400;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;
for m=1:runs
    m
    for i=1:N
        while true
            pop(i,:)=popMin+(popMax-popMin).*rand(1,D); 
            V(i,:)=vMin+(vMax-vMin).*rand(1,D);
            x=pop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            if x(D+1)>=DMin && x(D+1)<=DMax, break;end
        end
        fitness(i)=power_allocation_3D(x);
    end
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
            fitness(j)=power_allocation_3D(x);
            if fitness(j)<fitnesspbest(j)
                pbest(j,:)=pop(j,:);
                fitnesspbest(j)=fitness(j);
            end

            if fitness(j)<fitnessgbest
                gbest=pop(j,:);
                fitnessgbest=fitness(j);
            end

        end
%         w=w-dw;
        yy(i+1)=fitnessgbest;  
    end
%     plot(yy,'r');
%     hold on;
    COST(m,1)=fitnessgbest;
    COST(m,1);
    COSTS(m,:)= yy;
%     xlabel('Generation');
%     ylabel('fitness');
%     title('Power Allocation PSO-3D ');
end
% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('PSO-3D MEAN');
save('PSO_3D','COSTS')

