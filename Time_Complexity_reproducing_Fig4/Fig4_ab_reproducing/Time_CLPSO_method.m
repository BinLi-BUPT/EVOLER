% Copyright 2022, All Rights Reserved
% Code by Bin Li, Ziping Wei

clc;
clear all;
close all
D=30;
ps=100;
num_all = 100;
popMax= 20;
popMin=-popMax;
maxGen=2600;
des_val= 1e-15;  %If the output result is less than the threshold 'des_val', the global optimum is found  
vMax=0.15*(popMax); 
vMin=(-1)*vMax;
PSO_1 = zeros(1,maxGen);
PSO_2 = zeros(1,maxGen);
PSO_3 = zeros(1,maxGen);
PSO_4 = zeros(1,maxGen);
PSO_5 = zeros(1,maxGen);
PSO_6 = zeros(1,maxGen);
PSO_7 = zeros(1,maxGen);
PSO_8 = zeros(1,maxGen);
PSO_9 = zeros(1,maxGen);
PSO_10 = zeros(1,maxGen);
PSO_11 = zeros(1,maxGen);
PSO_12 = zeros(1,maxGen);
downhill_13 = zeros(1,maxGen);
for run=1:num_all
    run
    M = (-2+4*rand(1,D));               % Shifted global optimum, range [-2,2]
    tic
    for i=1:ps
        initPop(i,:)=popMin+(popMax-popMin)*rand(1,D);
        initV(i,:)=vMin+(vMax-vMin)*rand(1,D);
        initFitness(i)=shifted_levy(initPop(i,:),M);%*****************************************************************
    end
    initpBest=initPop;
    initpBestFit=initFitness;
    [initgBestFit,index]=min(initpBestFit);
    initgBest=initpBest(index);
    time_ini = toc;
    

    %% CLPSO 5
    tic
    c=1.49445;
    w=0.9;
    dw=(0.9-0.4)/maxGen;
    maxFlag=7;
    pop=initPop;
    v=initV;
    fitness=initFitness;
    pBest=pop;
    pBestFit=initFitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index);
    yy5(1)=gBestFit;
    PDF_cal5 = 0;
    for i=1:ps
        p(i)=0.05+0.45*(exp(10*(i-1)/(ps-1))-1)/(exp(10)-1);
    end
    for i=1:ps
        f=i*ones(1,D);
        fi1=ceil(ps.*rand(1,D));
        fi2=ceil(ps.*rand(1,D));
        fi=(pBestFit(fi1)<pBestFit(fi2)).*fi1+(pBestFit(fi1)>=pBestFit(fi2)).*fi2;
        RAND=rand(1,D);
        IND=find(RAND<p(i));
        if length(RAND)>0, f(IND)=fi(IND); 
        else
            randD=ceil(D*rand);
            f(randD)=fi(randD);
        end
        target(i,:)=f;
        flag(i)=0;
    end   
    for k=1:maxGen
        for i=1:ps
            if flag(i)==maxFlag
                f=i*ones(1,D);
                fi1=ceil(ps.*rand(1,D));
                fi2=ceil(ps.*rand(1,D));
                fi=(pBestFit(fi1)<pBestFit(fi2)).*fi1+(pBestFit(fi1)>=pBestFit(fi2)).*fi2;
                RAND=rand(1,D);
                IND=find(RAND<p(i));
                if length(RAND)>0, f(IND)=fi(IND);
                else
                    randD=ceil(D*rand);
                    f(randD)=fi(randD);
                end   
                target(i,:)=f;
                flag(i)=0;
            end
            %%
            for d=1:D
                index=target(i,d);
                v(i,d)=w*v(i,d)+c*rand*(pBest(index,d)-pop(i,d));
                if v(i,d)>vMax, v(i,d)=vMax;end
                if v(i,d)<vMin, v(i,d)=vMin;end
                pop(i,d)=pop(i,d)+v(i,d);
                if pop(i,d)>popMax, pop(i,d)=popMax; end
                if pop(i,d)<popMin, pop(i,d)=popMin; end
            end
            fitness(i)=shifted_levy(pop(i,:),M);%****************************************************************
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
        yy5(k+1)=gBestFit;
        if gBestFit > des_val
            PDF_cal5 = PDF_cal5 + 1;
        end     
    end 
    PDF_num_CLPSO(run) = PDF_cal5;
    time5(run) = toc;

    if yy5(maxGen - 5) < des_val
        PSO_5(run) = 1;
    end
    ff_end5(run,:) = yy5;
end
ff_end5 = ff_end5(:,1:2600);

save CLPSO_data ff_end5

