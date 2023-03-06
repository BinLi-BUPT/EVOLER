clear all;
clc;
format  long;
runs=50;
maxGen=500;
ps=20;
D=5;%
popSum=1260;
popMax=[500,200,300,150,200];
popMin=[100,50,80,50,50];
DMin=50; DMax=120;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;

for m=1:runs
    m
   for i=1:ps
       while(1)
            initPop(i,:)=popMin+(popMax-popMin).*rand(1,D);
            x=initPop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            if x(D+1)<=DMax && x(D+1)>=DMin, break;end
       end
        initFitness(i)=power_allocation_6D(x);%*******************************
        initV(i,:)=vMin+(vMax-vMin).*rand(1,D);
    end
    initpBest=initPop;
    initpBest=initFitness;
    [initgBestFit,gBestIndex]=min(initFitness);
    initgBest=initPop(gBestIndex);
    %DPSO
    c1=2;c2=2;
    w=0.4;
    cv=0;
    cl=0.001;
    pop=initPop;
    v=initV;
    fitness=initFitness;
    pBest=initPop;
    pBestFit=initFitness;
    gBest=initgBest;
    gBestFit=initgBestFit;
    count=1;
    yy2(1)=gBestFit;
    for i=1:maxGen
        for j=1:ps
            v(j,:)=w*v(j,:)+c1*rand*(pBest(j,:)-pop(j,:))+c2*rand*(gBest-pop(j,:));
            for d=1:D
                if v(j,d)>vMax(d),  v(j,d)=vMax(d); end
                if v(j,d)<vMin(d),  v(j,d)=vMin(d); end
            end
            for d=1:D
                if rand<cv
                    v(j,d)=rand*vMax;
                end
            end
            pop(j,:)=pop(j,:)+v(j,:);
            for d=1:D
                if pop(j,d)>popMax(d), pop(j,d)=popMax(d);  end
                if pop(j,d)<popMin(d), pop(j,d)=popMin(d);  end
            end
            for d=1:D
                if rand<cl
                    pop(j,d)=popMin(d)+(popMax(d)-popMin(d))*rand;
                end
            end
            x=pop(j,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            fitness(j)=power_allocation_6D(x);%*******************************
            if fitness(j)<pBestFit(j)
                pBest(j,:)=pop(j,:);
                pBestFit(j)=fitness(j);
            end
            if fitness(j)<gBestFit
                gBest=pop(j,:);
                gBestFit=fitness(j);
            end
%             yy2(count)=gBestFit;
%             count=count+1;
        end
        yy2(i+1)=gBestFit;
%          w=w-dw;
    end
%     plot(yy2,'r');
%     hold on;
%     title('Power Allocation DPSO-6D ');
    COST(m,1)=gBestFit;
    COST(m,1)
    COSTS(m,:)=yy2;
end

% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('DPSO-6D MEAN');
save('DPSO_6D','COSTS')
