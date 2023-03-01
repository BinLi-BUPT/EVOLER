%CPSO-H-6D

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
   % CPSO-H
    m
    w=0.9;
    c1=1.49;c2=1.49;
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
        initFitness(i)=power_allocation_6D(x);
        initV(i,:)=vMin+(vMax-vMin).*rand(1,D);
    end
    pop=initPop;
    v=initV;
    fitness=initFitness;
    pBest=pop;
    pbest1=pBest;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index,:);
    gBest1=gBest;
    s_c1=1.49;s_c2=1.49;
    s_ps=sizepop/2;
    s_w=0.9;
    for i=1:s_ps
        while(1)
            s_initPop(i,:)=popMin+(popMax-popMin).*rand(1,D);
            x=s_initPop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            if x(D+1)<=DMax && x(D+1)>=DMin, break; end
        end
        s_initFitness(i)=power_allocation_6D(x);
        s_initV(i,:)=vMin+(vMax-vMin).*rand(1,D);
    end
    s_pop=s_initPop;
    s_v=s_initV;
    s_fitness=s_initFitness;
    s_pBest=s_pop;
    s_pBestFit=s_fitness;
    [s_gBestFit,index]=min(s_pBestFit);
    s_gBest=s_pBest(index,:);
    %% 
    yy1(1)=min(gBestFit,s_gBestFit);
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
                fit2=power_allocation_6D(y);
                if fit1<fit2
                    gBestFit=fit1;
                    gBest(j)=pBest(i,j);
                end 
            end
        end
        w=w-dw;
        bool=1;
        for iii=1:20
            randIndex=randi(s_ps/2);%?
            for j=1:D
                if s_pop(randIndex,j)~=s_gBest(j)
                    bool=0;
                    break;
                end
            end
        end
        s_pop(randIndex,:)=gBest;%
        for i=1:ps
            s_v(i,:)=s_w*s_v(i,:)+s_c1*rand*(s_pBest(i,:)-s_pop(i,:))+s_c2*rand*(s_gBest-s_pop(i,:));
            for d=1:D
                if s_v(i,d)>vMax(d), s_v(i,d)= vMax(d); end 
                if s_v(i,d)<vMin(d), s_v(i,d)= vMin(d);  end
            end
            s_pop(i,:)=s_pop(i,:)+s_v(i,:);
            for d=1:D
                if s_pop(i,d)>popMax(d),    s_pop(i,d)=popMax(d);   end
                if s_pop(i,d)<popMin(d),    s_pop(i,d)=popMin(d);   end
            end
            x=s_pop(i,:);
            x(D+1)=popSum;
            for d=1:D
                x(D+1)=x(D+1)-x(d);
            end
            s_fitness(i)=power_allocation_6D(x);%**************************
            if s_fitness(i)<s_pBestFit(i)
                s_pBestFit(i)=s_fitness(i);
                s_pBest(i,:)=s_pop(i,:);
            end
            if s_fitness(i)<s_gBestFit
                s_gBestFit=s_fitness(i);
                s_gBest=s_pop(i,:);
            end
        end
        s_w=s_w-dw;
        for d=1:D
            bool=1;
            for iii=1:20
                randIndex=randi(ps/2);
                if pop(randIndex,d)~=gBest(d)
                    bool=0;
                    break;
                end
            end
            pop(randIndex,d)=s_gBest(d);
        end
        %%
        yy1(k+1)=min(gBestFit,s_gBestFit);
    end
    COST(m,1)=min(gBestFit,s_gBestFit);
    COST(m,1)
    COSTS(m,:)=yy1;
%     plot(yy1,'r');
%     hold on;
%     xlabel('Generation');
%     ylabel('fitness');
%     title('Power Allocation CPSO-H-6D');
end
% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('CPSO-H-6D MEAN');
save('CPSO_H_6D','COSTS')
