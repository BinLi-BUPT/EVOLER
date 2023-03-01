%CPSO-S-6D
clc;
clear all;
format long;
runs=50;
D=5;
sizepop=20;
popSum=1260;
popMax=[500,200,300,150,200];
popMin=[100,50,80,50,50];
DMin=50; DMax=120;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;
maxGen=500;
c1=1.49;c2=1.49;
dw=(0.9-0.4)/maxGen;

for m=1:runs
    m
    ps=ceil(sizepop/D);
    w=0.9;
    for i=1:ps
        while(1)
            initPop(i,:)=popMin+(popMax-popMin).*rand(1,D);
            x=initPop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            if x(D+1)<=DMax && x(D+1)>=DMin, break; end
        end
        initFitness(i)=power_allocation_6D(x);%%%
        initV(i,:)=vMin+(vMax-vMin).*rand(1,D);
    end
   % CPSO-S
    c1=1.49;c2=1.49;
    w=0.9;
    pop=initPop;
    v=initV;
    fitness=initFitness;

    pBest=pop;
    pbest1=pBest;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index,:);
    gBest1=gBest;
    yy1(1)=gBestFit;
    count=1;
    for k=1:maxGen
        for j=1:D
            for i=1:ps
                v(i,j)=w*v(i,j)+c1*rand*(pBest(i,j)-pop(i,j))+c2*rand*(gBest(j)-pop(i,j));
                if v(i,j)>vMax(j), v(i,j)=vMax(j);end
                if v(i,j)<vMin(j), v(i,j)=vMin(j);end
                pop(i,j)=pop(i,j)+v(i,j);
                if pop(i,j)>popMax(j), pop(i,j)=popMax(j);end
                if pop(i,j)<popMin(j), pop(i,j)=popMin(j);end
                temp1=gBest;pBest1=gBest;
                temp1(j)=pop(i,j);pBest1(j)=pBest(j);
                x=temp1;    y=pBest1;
                x(D+1)=popSum;      y(D+1)=popSum; 
                for d=1:D
                    x(D+1)=x(D+1)-x(d); 
                    y(D+1)=y(D+1)-y(d); 
                end
                fit1=power_allocation_6D(x);
                fit2=power_allocation_6D(y);
                if fit1<fit2
                    pBestFit(i)=fit1;
                    pBest(i,j)=pop(i,j);
                end
                temp2=gBest;gBest2=gBest;
                temp2(j)=pBest(i,j);
                x=temp2;    y=gBest2;
                x(D+1)=popSum;      y(D+1)=popSum; 
                for d=1:D
                    x(D+1)=x(D+1)-x(d); 
                    y(D+1)=y(D+1)-y(d); 
                end
                fit1=power_allocation_6D(x);
                fit2=power_allocation_6D(y);%%%
                if fit1<fit2
                    gBestFit=fit1;
                    gBest(j)=pBest(i,j);
                end 
            end
        end
        yy1(count+1)=gBestFit;
        count=count+1;
        w=w-dw;
    end 
    COST(m,1)=gBestFit;
    COST(m,1)
    COSTS(m,:)=yy1;
%     plot(yy1,'r');
%     hold on 
%     xlabel('Generation');
%     ylabel('fitness');
%     title('Power Allocation CPSO-S-6D ');
end
% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('CPSO-S-6D MEAN'); 
save('CPSO_S_6D','COSTS')
        


