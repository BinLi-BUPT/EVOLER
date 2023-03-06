
clear all;
clc;
format long;
runs=50;
c1=2.05;c2=2.05;c3=2.0;
groupNum=4;
ps=5;
maxGen=1000;

dw=(0.9-0.4)/maxGen;
D=2;
popSum=850;
popMax=[600,200];
popMin=[100,50];
DMin=100;
DMax=400;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;
for m=1:runs
    m
    w=0.9;
    for i=1:groupNum
        for j=1:ps
            while(1)
                Pop(i,j,:)=popMin+(popMax-popMin).*rand(1,D);
                for d=1:D
                    if Pop(i,j,d)>popMax(d),    Pop(i,j,d)=popMax(d);   end
                    if Pop(i,j,d)<popMin(d),    Pop(i,j,d)=popMin(d);   end
                end
                V(i,j,:)=vMin+(vMax-vMin).*rand(1,D);
                for d=1:D
                    if V(i,j,d)>vMax(d),    V(i,j,d)=vMax(d);   end
                    if V(i,j,d)<vMin(d),    V(i,j,d)=vMin(d);   end
                end
                x=Pop(i,j,:);
                x(D+1)=popSum;
                for d=1:D
                    x(D+1)=x(D+1)-x(d);
                end
                if x(D+1)>=DMin && x(D+1)<=DMax, break;end 
            end
            Fitness(i,j)=power_allocation_3D(x);%*******************************
        end
        [GBestFit(i),gBestIndex(i)]=min(Fitness(i,:));
        GBest(i,:)=Pop(i,gBestIndex(i),:);
        PBest=Pop;
        PBestFit=Fitness;
    end
    [globalBestFit(1),index]=min(GBestFit(i));
    globalBest=GBestFit(index,:);


    for i=1:maxGen
        for j=1:groupNum-1
            for k=1:ps
                for d=1:D
                    V(j,k,d)=w*V(j,k,d)+c1.*rand.*(PBest(j,k,d)-Pop(j,k,d))+c2.*rand.*(GBest(j,d)-Pop(j,k,d));
                    if V(j,k,d)>vMax(d),    V(j,k,d)=vMax(d);   end
                    if V(j,k,d)<vMin(d),    V(j,k,d)=vMin(d);   end
                    Pop(j,k,d)=Pop(j,k,d)+V(j,k,d);
                    if Pop(j,k,d)>popMax(d),    Pop(j,k,d)=popMax(d);   end
                    if Pop(j,k,d)<popMin(d),    Pop(j,k,d)=popMin(d);   end
                end          
                x=Pop(j,k,:);
                x(D+1)=popSum;
                for d=1:D
                    x(D+1)=x(D+1)-x(d);
                end
                Fitness(j,k)=power_allocation_3D(x);%*****************************

                if Fitness(j,k)<PBestFit(j,k)
                    PBestFit(j,k)=Fitness(j,k);
                    PBest(j,k,:)=Pop(j,k,:);
                end
                if Fitness(j,k)<GBestFit(j)
                    GBestFit(j)=Fitness(j,k);
                    GBest(j,:)=Pop(j,k,:);
                end
            end
            fitnesses(j,i)=GBestFit(j);

        end
        [slaveBestFit,bestGroup]=min(GBestFit);
        slaveBest=GBest(bestGroup,:);
        if slaveBestFit>GBestFit(groupNum)
            phi1=1;phi2=0;
        end
        if slaveBestFit==GBestFit(groupNum)
            phi1=0.5;phi2=0.5;
        end
        if slaveBestFit<GBestFit(groupNum)
            phi1=0;phi2=1;
        end
        for k=1:ps
            for d=1:D
                V(groupNum,k,d)=w*V(groupNum,k,d)+c1*rand*(PBest(groupNum,k,d)-Pop(groupNum,k,d))+c2*rand*phi1*(GBest(groupNum,d)-Pop(groupNum,k,d))+c3*rand*phi2*(slaveBest(d)-Pop(groupNum,k,d));
                if V(groupNum,k,d)>vMax(d),  V(groupNum,k,d)=vMax(d);    end
                if V(groupNum,k,d)<vMin(d),  V(groupNum,k,d)=vMin(d);    end
                Pop(groupNum,k,d)=Pop(groupNum,k,d)+V(groupNum,k,d);
                if Pop(groupNum,k,d)>popMax(d), Pop(groupNum,k,d)=popMax(d);    end
                if Pop(groupNum,k,d)<popMin(d), Pop(groupNum,k,d)=popMin(d);    end
            end
            x=Pop(groupNum,k,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            fit=power_allocation_3D(x);%*****************************
            if fit<PBestFit(groupNum,k)
                PBestFit(groupNum,k)=fit;
                PBest(groupNum,k,:)=Pop(groupNum,k,:);
            end
            if fit < GBestFit(groupNum)
                GBestFit(groupNum)=fit;
                GBest(groupNum,:)=Pop(groupNum,k,:);
            end
        end
        fitnesses(groupNum,i)=GBestFit(groupNum);
        [globalBestFit(i+1),index]=min(fitnesses(:,i));
        w=w-dw;
    end
    COST(m,1)=min(fitnesses(:,i));
    COST(m,1)
    COSTS(m,:)=globalBestFit;
    plot(globalBestFit,'r');
%     hold on;
%     xlabel('Generation');
%     ylabel('fitness')
%     title('Power Allocation MPSO-3D ')
end
% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('MPSO-3D MEAN');
save('MCPSO_3D','COSTS')
