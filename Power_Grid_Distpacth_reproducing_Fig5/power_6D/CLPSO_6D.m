% CLPSO-6D
clc;
clear all;
format long;

runs=50;
D=5;
ps=20;
popSum=1260;
popMax=[500,200,300,150,200];
popMin=[100,50,80,50,50];
DMin=50; DMax=120;
vMax=0.2*(popMax-popMin); vMin=(-1)*vMax;
maxGen=500;
c=1.49445;
dw=(0.9-0.4)/maxGen;
maxFlag=7;
for m=1:runs
    m
    w=0.9;
    for i=1:ps
        while true
            initPop(i,:)=popMin+(popMax-popMin).*rand(1,D);
            x=initPop(i,:);
            x(D+1)=popSum;
            for d=1:D,  x(D+1)=x(D+1)-x(d);     end
            if x(D+1)>=DMin && x(D+1)<=DMax, break; end
        end
        initFitness(i)=power_allocation_6D(x);%*******
        initV(i,:)=vMin+(vMax-vMin).*rand(1,D);
    end
    pop=initPop;
    v=initV;
    fitness=initFitness;
    pBest=pop;
    pBestFit=initFitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index);
    yy1(1)=gBestFit;
    for i=1:ps
        p(i)=0.05+0.45*(exp(10*(i-1)/(ps-1))-1)/(exp(10)-1);
    end
    for i=1:ps
        f=i*ones(1,D);
        fi1=ceil(ps.*rand(1,D));
        fi2=ceil(ps.*rand(1,D));
        fi=(pBestFit(fi1)<pBestFit(fi2)).*fi1+(pBestFit(fi1)>=pBestFit(fi2)).*fi2;
        for d=1:D
            if rand<=p(i),f(d)=fi(d);end
        end
        if f==i.*ones(1,D)
            randD=ceil(D*rand);
            f(randD)=fi(randD);
        end
        target(i,:)=f;
        flag(i)=0;
    end
    %%

    for k=1:maxGen
        for i=1:ps
            if flag(i)==maxFlag
                f=i*ones(1,D);
                fi1=ceil(ps.*rand(1,D));
                fi2=ceil(ps.*rand(1,D));
                fi=(pBestFit(fi1)<pBestFit(fi2)).*fi1+(pBestFit(fi1)>=pBestFit(fi2)).*fi2;
                for d=1:D
                    if rand<=p(i),f(d)=fi(d);end
                end
                if f==i.*ones(1,D)
                    randD=ceil(D*rand);
                    f(randD)=fi(randD);
                end
                target(i,:)=f;
                flag(i)=0;
            end
            for d=1:D
                index=target(i,d);
                v(i,d)=w*v(i,d)+c*rand*(pBest(index,d)-pop(i,d));
                if v(i,d)>vMax(d), v(i,d)=vMax(d);end
                if v(i,d)<vMin(d), v(i,d)=vMin(d);end
                pop(i,d)=pop(i,d)+v(i,d);
                if pop(i,d)>popMax(d), pop(i,d)=popMax(d); end
                if pop(i,d)<popMin(d), pop(i,d)=popMin(d); end
            end
            x=pop(i,:);
            x(D+1)=popSum;
            for d=1:D,  x(D+1)=x(D+1)-x(d);     end
            fitness(i)=power_allocation_6D(x(1,:));%*******
            if fitness(i)<pBestFit(i)
                pBest(i,:)=pop(i,:);
                pBestFit(i)=fitness(i);
                flag(i)=0;
            else, flag(i)=flag(i)+1;
            end
            if pBestFit(i)<gBestFit
                gBestFit=pBestFit(i);
                gBest=pBest(i,:);
            end
        end
        w=w-dw;
        yy1(k+1)=gBestFit;
    end
    COST(m,:)=gBestFit;
    COST(m)
    COSTS(m,:)=yy1;
%     plot(yy1,'r');
%     hold on;
%     xlabel('Generation');
%     ylabel('fitness');
%     title('Power Allocation CLPSO-6D ');
end
% hold off;
% figure;
% for i=1:(maxGen+1)
%     AVE(i)=mean(COSTS(:,i));
% end
% plot(AVE,'r');
% title('CLPSO-6D MEAN');
save('CLPSO_6D_50times','COSTS')
