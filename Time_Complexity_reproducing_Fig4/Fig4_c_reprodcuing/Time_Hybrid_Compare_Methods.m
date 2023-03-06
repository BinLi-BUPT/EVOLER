% Copyright 2022, All Rights Reserved
% Code by Bin Li, Ziping Wei
% 
% % Several Comparing methods
% Standard PSO、CLPSO、CLPSO-H、CLPSO-S、ML-PSO、DPSO、NPSO、GCPSO、MPCPSO、DNSPSO、Local PSO、Downhill
clc;
clear all;
close all
%% 
D=30;
ps=100;
num_all = 20;
maxGen=900;
popMax= 5;
popMin=-popMax;
des_val= 8e-4;  %If the output result is less than the threshold 'des_val', the global optimum is found  
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
    tic
    for i=1:ps
        initPop(i,:)=popMin+(popMax-popMin)*rand(1,D);
        initV(i,:)=vMin+(vMax-vMin)*rand(1,D);
        initFitness(i)=hybrid_func1(initPop(i,:));%*****************************************************************
    end
    initpBest=initPop;
    initpBestFit=initFitness;
    [initgBestFit,index]=min(initpBestFit);
    initgBest=initpBest(index);
    time_ini = toc;
    
    %% Standard PSO 1
    tic
    c1=1.4962;c2=1.4962;
    w=0.7968;
    pop=initPop;
    v=initV;
    fitness=initFitness;
    pBest=pop;
    pBestFit=initFitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index);
    yy1(1)=gBestFit;
    PDF_cal = 0;
    for i=1:maxGen
        for j=1:ps
            v(j,:)=w*v(j,:)+c1*rand(1,D).*(pBest(j,:)-pop(j,:))+c2*rand(1,D).*(gBest-pop(j,:));
            v(j,:)=w*v(j,:)+c1*rand.*(pBest(j,:)-pop(j,:))+c2*rand.*(gBest-pop(j,:));
            v(j,find(v(j,:)>vMax))=vMax;
            v(j,find(v(j,:)<vMin))=vMin;
            pop(j,:)=pop(j,:)+v(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            fitness(j)=hybrid_func1(pop(j,:));%************************************************************************
            if fitness(j)<pBestFit(j)
                pBest(j,:)=pop(j,:);
                pBestFit(j)=fitness(j);
            end
            if fitness(j)<gBestFit
                gBest=pop(j,:);
                gBestFit=fitness(j);
            end
        end
        yy1(i+1)=gBestFit;
        if gBestFit > des_val
            PDF_cal = PDF_cal + 1;
        end     
    end
    time1(run) = toc;
    
    
    PDF_num_SPSO(run) = PDF_cal;
    if yy1(maxGen - 5) < des_val
        PSO_1(run) = 1;
    end
    ff_end1(run,:) = yy1;
    %% DPSO 2
    tic
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
    yy2(1)=gBestFit;
    PDF_cal2 = 0;
    for i=1:maxGen
        for j=1:ps
            v(j,:)=w*v(j,:)+c1*rand(1,D).*(pBest(j,:)-pop(j,:))+c2*rand(1,D).*(gBest-pop(j,:));
            v(j,find(v(j,:)>vMax))=vMax;
            v(j,find(v(j,:)<vMin))=vMin;
            RAND=rand(1,D);
            IND=find(RAND<cv);
            if length(IND)>0, v(j,IND)=vMax.*rand(1,length(IND)); end
            pop(j,:)=pop(j,:)+v(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            RAND=rand(1,D);
            IND=find(RAND<cl);
            if length(IND)>0, pop(j,IND)=popMin+(popMax-popMin).*rand(1,length(IND)); end 
            fitness(j)=hybrid_func1(pop(j,:));%**********************************************************************
            if fitness(j)<pBestFit(j)
                pBest(j,:)=pop(j,:);
                pBestFit(j)=fitness(j);
            end
            if fitness(j)<gBestFit
                gBest=pop(j,:);
                gBestFit=fitness(j);
            end
        end
        yy2(i+1)=gBestFit;
        if gBestFit > des_val
            PDF_cal2 = PDF_cal2 + 1;
        end     
    end
    time2(run) = toc;

    PDF_num_DPSO(run) = PDF_cal2;
    if yy2(maxGen - 5) < des_val
        PSO_2(run) = 1;
    end
    ff_end2(run,:) = yy2;
    %% CPSO-S 3
    tic
    c1=1.49;c2=1.49;
    w=0.9;
    dw=(0.9-0.4)/maxGen;
    pop=initPop(1:ceil(ps/D),:);
    v=initV(1:ceil(ps/D),:);
    fitness=initFitness(1:ceil(ps/D));
    pBest=pop;
    pbest1=pBest;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index,:);
    gBest1=gBest;
    yy3(1)=gBestFit;
    count=1;
    PDF_cal3 = 0;
    for k=1:maxGen
        for j=1:D
            for i=1:ceil(ps/D)
                v(i,j)=w*v(i,j)+c1*rand*(pBest(i,j)-pop(i,j))+c2*rand*(gBest(j)-pop(i,j));
                if v(i,j)>vMax, v(i,j)=vMax;end
                if v(i,j)<vMin, v(i,j)=vMin;end
                pop(i,j)=pop(i,j)+v(i,j);
                if pop(i,j)>popMax, pop(i,j)=popMax;end
                if pop(i,j)<popMin, pop(i,j)=popMin;end
                temp1=gBest;pBest1=gBest;
                temp1(j)=pop(i,j);pBest1(j)=pBest(j);
                 fit1=hybrid_func1(temp1);
                 fit2=hybrid_func1(pBest1);%***************************************************************
                if fit1<fit2
                    pBestFit(i)=fit1;
                    pBest(i,j)=pop(i,j);
                end
                
                if pBestFit(i)<gBestFit
                    gBestFit=pBestFit(i);
                    gBest(j)=pBest(i,j);
                end 
            end
        end
        yy3(k+1)=gBestFit;
        count=count+1;
        w=w-dw;
        if gBestFit > des_val
            PDF_cal3 = PDF_cal3 + 1;
        end     
    end
    PDF_num_CPSO(run) = PDF_cal3;
    if yy3(maxGen - 5) < des_val
        PSO_3(run) = 1;
    end
    ff_end3(run,:) = yy3;
    t3=clock;
    t3=t3(4:6);
    time3(run) = toc;

     %% CPSO-H 4
    tic
    c1=1.49;c2=1.49;
    w=0.9;
    dw=(0.9-0.4)/maxGen;
    p_s=ceil(ps/D);
    pop=initPop(1:p_s,:);
    v=initV(1:p_s,:);
    fitness=initFitness(1:p_s);
    pBest=pop;
    pbest1=pBest;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index,:);
    gBest1=gBest;
    
    s_c1=1.49;s_c2=1.49;
    s_ps=ps/2;
    s_w=0.9;
    s_dw=(0.9-0.4)/maxGen;
    s_pop=initPop((p_s+1):ps,:);
    s_v=initV((p_s+1):ps,:);
    s_fitness=initFitness((p_s+1):ps);
    s_pBest=s_pop;
    s_pBestFit=s_fitness;
    [s_gBestFit,index]=min(s_pBestFit);
    s_gBest=s_pBest(index,:);
    yy4(1)=min(gBestFit,s_gBestFit);
    PDF_cal4 = 0;
    for k=1:maxGen
        for j=1:D
            for i=1:p_s
                v(i,j)=w*v(i,j)+c1*rand*(pBest(i,j)-pop(i,j))+c2*rand*(gBest(j)-pop(i,j));
                if v(i,j)>vMax, v(i,j)=vMax;end
                if v(i,j)<vMin, v(i,j)=vMin;end
                pop(i,j)=pop(i,j)+v(i,j);
                if pop(i,j)>popMax, pop(i,j)=popMax;end
                if pop(i,j)<popMin, pop(i,j)=popMin;end
                temp1=gBest;pBest1=gBest;
                temp1(j)=pop(i,j);pBest1(j)=pBest(j);
                fit1=hybrid_func1(temp1);fit2=hybrid_func1(pBest1);%****************************************************
                if fit1<fit2
                    pBestFit(i)=fit1;
                    pBest(i,j)=pop(i,j);
                end
                
                if pBestFit(i)<gBestFit
                    gBestFit=pBestFit(i);
                    gBest(j)=pBest(i,j);
                end 
            end
        end
        w=w-dw;
        randIndex=randi(ceil(s_ps/2));
        if s_pop(randIndex,:)==s_gBest
            if randIndex==1, randIndex=randIndex+1;
            else, randIndex=randIndex-1;
            end
        end
 
        s_pop(randIndex,:)=gBest;%
        for i=1:p_s
            s_v(i,:)=s_w*s_v(i,:)+s_c1*rand(1,D).*(s_pBest(i,:)-s_pop(i,:))+s_c2*rand(1,D).*(s_gBest-s_pop(i,:));
            s_v(i,find(s_v(i,:)>vMax))=vMax;
            s_v(i,find(s_v(i,:)<vMin))=vMin;
            s_pop(i,:)=s_pop(i,:)+s_v(i,:);
            s_pop(i,find(s_pop(i,:)>popMax))=popMax;
            s_pop(i,find(s_pop(i,:)<popMin))=popMin;
            s_fitness(i)=hybrid_func1(s_pop(i,:));%*************************************************************
            if s_fitness(i)<s_pBestFit(i)
                s_pBestFit(i)=s_fitness(i);
                s_pBest(i,:)=s_pop(i,:);
            end
            if s_fitness(i)<s_gBestFit
                s_gBestFit=s_fitness(i);
                s_gBest=s_pop(i,:);
            end
        end
        s_w=s_w-s_dw;
        for d=1:D
            randIndex=randi(ceil(p_s/2));
            if pop(randIndex,d)==gBest(d) 
                if randIndex==1, randIndex=randIndex+1;
                else, randIndex=randIndex-1;
                end
            end
            pop(randIndex,d)=s_gBest(d);
        end
        %
        yy4(k+1)=min(gBestFit,s_gBestFit);
        if min(gBestFit,s_gBestFit) > des_val
            PDF_cal4 = PDF_cal4 + 1;
        end     
    end  %% 
    time4(run) = toc;

    PDF_num_CPSOH(run) = PDF_cal4;
    if yy4(maxGen - 5) < des_val
        PSO_4(run) = 1;        
    end
    ff_end4(run,:) = yy4;
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
            fitness(i)=hybrid_func1(pop(i,:));%****************************************************************
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
    tic
    %% MCPSO 6
    c1=2.05;c2=2.05;c3=2.0;
    groupNum=5;%
    p_s=ps/groupNum;
    w=0.9;
    dw=(0.9-0.4)/maxGen;
    PDF_cal6 = 0;
    V = [];
    Pop = [];
    Fitness = [];
    for i=1:groupNum
        Pop(i,:,:)=initPop((i-1)*p_s+(1:p_s),:);
        V(i,:,:)=initV((i-1)*p_s+(1:p_s),:);
        Fitness(i,:)=initFitness((i-1)*p_s+(1:p_s));
        [GBestFit(i),gBestIndex(i)]=min(Fitness(i,:));
        GBest(i,:)=Pop(i,gBestIndex(i),:);
        PBest=Pop;
        PBestFit=Fitness;
    end
    [globalBestFit,index]=min(GBestFit);%%%%%%%%%%%%%%
    globalBest=GBest(index,:);%%%%%%%%%%%%
    yy6(1)=globalBestFit;
    for i=1:maxGen
        for j=1:groupNum-1
            for k=1:p_s
                for d=1:D
                    V(j,k,d)=w*V(j,k,d)+c1.*rand.*(PBest(j,k,d)-Pop(j,k,d))+c2.*rand.*(GBest(j,d)-Pop(j,k,d));
                end
                V(j,k,find(V(j,k,:)>vMax))=vMax;
                V(j,k,find(V(j,k,:)<vMin))=vMin;           
                Pop(j,k,:)=Pop(j,k,:)+V(j,k,:);
                Pop(j,k,find(Pop(j,k,:)>popMax))=popMax;
                Pop(j,k,find(Pop(j,k,:)<popMin))=popMin;
                Fitness(j,k)=hybrid_func1(reshape(Pop(j,k,:),1,D));%******************************************************************************
 
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
        [slaveBestFit,bestGroup]=min(GBestFit(1:groupNum-1));%%%%%%%%%%%%%%%%%%%
        slaveBest=GBest(bestGroup,:);%%%%%%%%%%%%%%%%%%%%%
 
        if slaveBestFit>GBestFit(groupNum)
            phi1=1;phi2=0;
        end
        if slaveBestFit==GBestFit(groupNum)
            phi1=0.5;phi2=0.5;
        end
        if slaveBestFit<GBestFit(groupNum)
            phi1=0;phi2=1;
        end
        for k=1:p_s
            for d=1:D
                V(groupNum,k,d)=w*V(groupNum,k,d)+c1*rand*(PBest(groupNum,k,d)-Pop(groupNum,k,d))+c2*rand*phi1*(GBest(groupNum,d)-Pop(groupNum,k,d))+c3*rand*phi2*(slaveBest(d)-Pop(groupNum,k,d));%%%%%%%%%%%%%%
            end
            V(groupNum,k,find(V(groupNum,k,:)>vMax))=vMax;
            V(groupNum,k,find(V(groupNum,k,:)<vMin))=vMin;
            Pop(groupNum,k,:)=Pop(groupNum,k,:)+V(groupNum,k,:);
            Pop(groupNum,k,find(Pop(groupNum,k,:)>popMax))=popMax;
            Pop(groupNum,k,Pop(groupNum,k,:)<popMin)=popMin;
            fit=hybrid_func1(reshape(Pop(groupNum,k,:),1,D));%*************************************************************************************
            Fitness(groupNum,k)=fit;
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
        [globalBestFit(i),index]=min(fitnesses(:,i));
        yy6(i+1)=min(fitnesses(:,i));
        w=w-dw;
        if min(fitnesses(:,i)) > des_val
            PDF_cal6 = PDF_cal6 + 1;
        end
    end  %%
    time6(run) = toc;

    PDF_num_MLPSO(run) = PDF_cal6;
    if yy6(maxGen - 5) < des_val
        PSO_6(run) = 1;
    end
    ff_end6(run,:) = yy6;
    %% NPSO 7
    tic
    c1=2;c2=2;
    pop=initPop;  
    V =initV;
    fitness=initFitness;
    gbest=initgBest;
    fitnessgbest=initgBestFit;
    [~, bestindex2]=max(initFitness);
    g_id=initPop(bestindex2,:);
    p_id=initPop;
    pbest=pop;
    fitnesspbest=fitness;
    PDF_cal7 = 0;
    yy7(1)=fitnessgbest;
    for i=1:maxGen
        for j=1:ps
            V(j,:)=V(j,:)+c1*rand(1,D).*(pop(j,:)-p_id(j,:))+c2*rand(1,D).*(pop(j,:)-g_id);
            V(j,find(V(j,:)>vMax))=vMax;
            V(j,find(V(j,:)<vMin))=vMin;
            pop(j,:) = pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            
            fitness(j)=hybrid_func1(pop(j,:));
            if fitness(j) < fitnesspbest(j)
                fitnesspbest(j)=fitness(j);
                pbest(j,:)=pop(j,:);
            else
                p_id(j,:) = pop(j,:);      
            end
            if fitness(j) < fitnessgbest
                fitnessgbest=fitness(j);
                gbest=pop(j,:);
            else
                g_id = pop(j,:);
            end     
        end
        yy7(i+1)=fitnessgbest;
        if fitnessgbest > des_val
            PDF_cal7 = PDF_cal7 + 1;
        end
    end  %%
        time7(run) = toc;

    PDF_num_NPSO(run) = PDF_cal7;
    if yy7(maxGen - 5) < des_val
        PSO_7(run) = 1;
    end
    ff_end7(run,:) = yy7;
     %% GCPSO 8
     tic
     npo = 2.1; fc = 5; sc = 15; c1=1.49; c2=1.49; w=0.72;
    fitness=initFitness;
    pop=initPop;
    gbest=initgBest;
    fitnessgbest=initgBestFit;
    V =initV;
    pbest=pop;
    fitnesspbest=fitness;
    gbest_all = zeros(1,maxGen);
    yy8(1) = initgBestFit;
    PDF_cal8 = 0;
    for i=1:maxGen
        for j=1:ps
            V(j,:)=w*V(j,:)+c1*rand(1,D).*(pbest(j,:)-pop(j,:))+c2*rand(1,D).*(gbest-pop(j,:));
            V(j,find(V(j,:)>vMax))=vMax;
            V(j,find(V(j,:)<vMin))=vMin;
            
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            
            fitness(j)=hybrid_func1(pop(j,:));%*****************************************************************
            if fitness(j)<fitnesspbest(j)
                pbest(j,:)=pop(j,:);
                fitnesspbest(j)=fitness(j);
            end
            if fitness(j)<fitnessgbest
                gbest=pop(j,:);
                fitnessgbest=fitness(j);
                gbest_all(i) = hybrid_func1(gbest);%*****************************************************************
                V_tao = - pop(j,:) + gbest + w*V(j,:)+ npo*(1-rand);
                pop(j,:) = gbest + w*V_tao + npo*(1-2*rand);
                gbest=pop(j,:);
            end
            if i > fc
                test_vec = length(unique(gbest_all(i-fc:i)));
                if test_vec == 1
                    npo =0.5*npo;
                end
                if i > sc
                    test_vec = length(unique(gbest_all(i-fc:i)));
                    if test_vec > 1
                        npo =2*npo;
                    end
                end
            end
        end
        yy8(i+1)=fitnessgbest;
        if fitnessgbest > des_val
            PDF_cal8 = PDF_cal8 + 1;
        end
    end  %
    time8(run) = toc;

    PDF_num_GCPSO(run) = PDF_cal8;
    if yy8(maxGen - 5) < des_val
        PSO_8(run) = 1;
    end
    ff_end8(run,:) = yy8;
    %% MPCPSO 10
    tic
    t=0.5;  MAXT=5; F=0.5;
    N1=t*ps;    N2=ps-N1;
    c=1.49445;  w=0.729;
    vmax=0.2*popMax; vmin=0.2*popMin;
    tournamentsize=5;
    pop=initPop;
    v=initV;
    fitness=initFitness;
    me1=sum(fitness);
    pbestfitness=fitness;
    pbest=pop;
    [sortfitness,sortindex]=sort(pbestfitness,2);
    gbestfitness=sortfitness(1);
    gbest=pbest(sortindex(1),:);
    T=0;
    yy10(1)=gbestfitness;
    PDF_cal10 = 0;
    for i=1:maxGen
        is_des=0;% 
         [sortfitness,sortindex]=sort(pbestfitness,2);
        EP=sortindex(1:N1);
        GP=sortindex(N1+1:end);
         segmentnum=ceil((i/maxGen)*D);
        ds=ceil(D/segmentnum);
        y=floor(D/segmentnum);
         sumfitness=0;
        for j=1:N1
            sumfitness=sumfitness+pbestfitness(GP(j));
        end
        phi=-1*me1*sumfitness/N1;
        phi=1/((1+exp(phi))^i);% 
        for j=1:N2
             index=GP(j);
             for dd=1:ds:D
                sum1=0;
                sum2=zeros(1,D);
                pl=zeros(1,D);
                for ii=1:tournamentsize
                    tournamentindex=EP(randi(N1)); 
                    sum1=sum1+pbest(tournamentindex);
                    sum2=sum2+pbestfitness(tournamentindex)*pbest(tournamentindex,:);
                end
                pl=sum2/(sum1*tournamentsize);
                r1=rand;r2=rand;
                if dd+ds-1>D
                    for d=(D-y+1):D
                        v(index,d)=w*v(index,d)+c*r1*(pl(d)-pop(index,d))+phi*r2*(gbest(d)-pop(index,d));
                        if v(index,d)<vmin, v(index,d)=vmin; end
                        if v(index,d)>vmax, v(index,d)=vmax; end
                        pop(index,d)=pop(index,d)+v(index,d);
                        if pop(index,d)<popMin, pop(index,d)=popMin; end
                        if pop(index,d)>popMax, pop(index,d)=popMax; end
                    end
                else
                    for d=dd:(dd+ds-1)
                        v(index,d)=w*v(index,d)+c*r1*(pl(d)-pop(index,d))+phi*r2*(gbest(d)-pop(index,d));
                        if v(index,d)<vmin, v(index,d)=vmin; end
                        if v(index,d)>vmax, v(index,d)=vmax; end
                        pop(index,d)=pop(index,d)+v(index,d);
                        if pop(index,d)<popMin, pop(index,d)=popMin; end
                        if pop(index,d)>popMax, pop(index,d)=popMax; end
                    end
                end
            end
            fitness(index)=hybrid_func1(pop(index,:));%*****************************************************************
            if fitness(index)<pbestfitness(index)
                pbestfitness(index)=fitness(index);
                pbest(index,:)=pop(index,:);
            end
            if pbestfitness(index)<gbestfitness
                gbestfitness=pbestfitness(index);
                gbest=pbest(index,:);
                is_des=1;T=0;
            end
        end
         
         if T<MAXT
             sumfitness=0;
            for j=1:N1
                sumfitness=sumfitness+pbestfitness(EP(j));
            end
            phi=-1*me1*sumfitness/N1;
            phi=1/((1+exp(phi))^i); 
             sumpop=zeros(1,D);
            for j=1:N1
                index=EP(j);
                sumpop=sumpop+pbest(index,:);
            end
            for j=1:N1
                index=EP(j);
                sumD=sum(pbest(index,:));
                 for d=1:D
                    lambda=1/(exp(sumD/D-pbest(index,d))+1);
                    r=rand;
                    mbest(d)=lambda*r/D*sumD+(1-lambda)*(1-r)*sumpop(d);
                    v(index,d)=w*v(index,d)+c*rand*(mbest(d)-pop(index,d));
                    if v(index,d)<vmin, v(index,d)=vmin; end
                    if v(index,d)>vmax, v(index,d)=vmax; end
                    pop(index,d)=phi*pop(index,d)+v(index,d);
                    if pop(index,d)<popMin, pop(index,d)=popMin; end
                    if pop(index,d)>popMax, pop(index,d)=popMax; end
                 end
                 fitness(index)=hybrid_func1(pop(index,:));%*****************************************************************
                if fitness(index)<pbestfitness(index)
                    pbestfitness(index)=fitness(index);
                    pbest(index,:)=pop(index,:);
                end
                if pbestfitness(index)<gbestfitness
                    gbestfitness=pbestfitness(index);
                    gbest=pbest(index,:);
                    is_des=1;T=0;
                end
            end
        else
            T=0;
            for j=1:N1
                index=EP(j);
                index1=randi(ps);
                index2=randi(ps);
                while(index1==index || index2==index || index1==index2)
                    index1=randi(ps);
                    index2=randi(ps);
                end
                pop(index,:)=gbest+F*(pbest(index1,:)-pbest(index2,:));
                pop(index,find(pop(index,:)>popMax))=popMax;
                pop(index,find(pop(index,:)<popMin))=popMin;
                 fitness(index)=hybrid_func1(pop(index,:));%*****************************************************************
                if fitness(index)<pbestfitness(index)
                    pbestfitness(index)=fitness(index);
                    pbest(index,:)=pop(index,:);
                end
                if pbestfitness(index)<gbestfitness
                    gbestfitness=pbestfitness(index);
                    gbest=pbest(index,:);
                    is_des=1;T=0;
                end
            end
         end
         if is_des==0, T=T+1; end
         result(i)=gbestfitness;
        yy10(i+1)=gbestfitness;
        if gbestfitness > des_val
            PDF_cal10 = PDF_cal10 + 1;
        end
    end
    time9(run) = toc;

    PDF_num_MPCPSO(run) = PDF_cal10;
     if yy10(maxGen - 5) < des_val
        PSO_10(run) = 1;
    end
    ff_end10(run,:) = yy10;
    %% DNSPSO 11
    tic
    PDF_cal11 = 0;
    tran_p=0.9;
    state=1;
    neighborsize=5;
    F=0.5;
    CR=0.9;
    dF=(1-0.5)/maxGen;
    dCR=(0.4-0.9)/maxGen;
    pop=zeros(ps,D);
    V=zeros(ps,D);
    fitness=zeros(1,ps);
    meandistance=zeros(1,ps);
    evfactor=zeros(1,ps);
    u=zeros(1,D);
    distance=zeros(1,ps);
    neighborindex=zeros(ps,neighborsize);
    status=zeros(1,ps);
    pbest2=zeros(1,D);
    gbest2=zeros(1,D);
    result=zeros(1,maxGen);
    pop=initPop;
    V=initV;
    fitness=initFitness;
    pbest=pop;
    pbestfitness=fitness;
    [gbestfitness,bestindex]=min(pbestfitness);
    gbest=pbest(bestindex,:);
    yy11(1)=gbestfitness;
    for i=1:maxGen
        i;
        for j=1:ps
            fitness(j)=hybrid_func1(pop(j,:));%***********************
            if fitness(j)<pbestfitness(j)
                pbestfitness(j)=fitness(j);
                pbest(j,:)=pop(j,:);
            end
            if pbestfitness(j)<gbestfitness
                gbestfitness=pbestfitness(j);
                gbest=pbest(j,:);
            end
        end
        for j=1:ps
            distance=0;
            for k=1:ps
                distance=distance+norm(pop(k,:)-pop(j,:));
            end
            meandistance(j)=distance/(ps-1);
        end
        distancemin=min(meandistance);
        distancemax=max(meandistance);
        for j=1:ps
            evfactor(j)=(meandistance(j)-distancemin)/(distancemax-distancemin);
            if evfactor(j)<0.25,status(j)=1;
            elseif evfactor(j)<0.5, status(j)=2;
            elseif evfactor(j)<0.75, status(j)=3;
            else,status(j)=4;
            end
        end
        for j=1:ps
            for d=1:D
                r11=randi(ps);
                r22=randi(ps);
                u(d)=pbest(j,d)+F*(pbest(r11,d)-pbest(r22,d));
                if u(d)<popMin,u(d)=popMin;end
                if u(d)>popMax,u(d)=popMax;end
                if rand>CR, u(d)=pbest(j,d);end
            end
            fitness_u=hybrid_func1(u);%**********************************
            if fitness_u<pbestfitness(j)
                pbest(j,:)=u;
                pbestfitness(j)=fitness_u;
            end
            if fitness_u<gbestfitness
                gbest=u;
                gbestfitness=fitness_u;
            end
        end
        result(i)=gbestfitness;
        for j=1:ps
            for k=1:ps
                distance(k)=norm(pbest(j,:)-pbest(k,:));
            end
            [distance2,sortindex]=sort(distance,2);
            neighborindex(j,:)=sortindex(1:neighborsize);
        end
        for k=1:ps
            distance(k)=norm(gbest-pbest(k,:));
        end
        [distance2,sortindex]=sort(distance,2);
        gneighborindex=sortindex(1:neighborsize);
        for j=1:ps
            r=rand;
            if status(j)==1
                if r>tran_p, status(j)=2;end
            elseif status(j)==2
                if r<(1-tran_p)/2, status(j)=1;
                elseif r>(1+tran_p)/2, status(j)=3;
                end
            elseif status(j)==3
                if r<(1-tran_p)/2, status(j)=2;
                elseif r>(1+tran_p)/2, status(j)=4;
                end
            else
                if r<(1-tran_p), status(j)=3;end
            end
        end
        w=0.5*evfactor+0.4;
        for j=1:ps
            if status(j)==1
                c11=2;c22=2;
                for d=1:D
                    pbest2(d)=pbest(j,d);
                end
                for d=1:D
                    gbest2(d)=gbest(d);
                end
            elseif status(j)==2
                c11=2.1;c22=1.9;
                r11=randi(ps);
                r22=randi(neighborsize);
                index=neighborindex(r11,r22);
                for d=1:D
                    pbest2(d)=pbest(index,d);
                end
                for d=1:D
                    gbest2(d)=gbest(d);
                end
            elseif status(j)==3
                c11=2.2;c22=1.8;
                for d=1:D
                    pbest2(d)=pbest(j,d);
                end
                r22=randi(neighborsize);
                index=gneighborindex(r22);
                for d=1:D
                    gbest2(d)=pbest(index,d);
                end
            else
                c11=1.8;c22=2.2;
                r11=randi(ps);
                r22=randi(neighborsize);
                index=neighborindex(r11,r22);
                for d=1:D
                    pbest2(d)=pbest(index,d);
                end
                r22=randi(neighborsize);
                index=gneighborindex(r22);
                for d=1:D
                    gbest2(d)=pbest(index,d);
                end
            end
            for d=1:D
                V(j,d)=w(j)*V(j,d)+c11*rand*(pbest2(d)-pop(j,d))+c22*rand*(gbest2(d)-pop(j,d));
                if V(j,d)<vMin,V(j,d)=vMin;end
                if V(j,d)>vMax,V(j,d)=vMax;end
                pop(j,d)=pop(j,d)+V(j,d);
                if pop(j,d)<popMin,pop(j,d)=popMin;end
                if pop(j,d)>popMax,pop(j,d)=popMax;end
            end
        end
        F=F+dF;
        CR=CR+dCR;
        yy11(i+1)=gbestfitness;
        if gbestfitness > des_val
            PDF_cal11 = PDF_cal11 + 1;
        end
    end
    time10(run) = toc;

    PDF_num_DNSPSO(run) = PDF_cal11;
    if yy11(maxGen - 5) < des_val
        PSO_11(run) = 1;
    end
    ff_end11(run,:) = yy11;
    t11=clock;
    t11=t11(4:6);
    %% Local PSO 12
    tic
    c1=1.4962; c2=1.4962; w=0.7968;
    pop=initPop;
    V=initV;
    fitness=initFitness;
    pBest=pop;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index);
    yy12(1)=gBestFit;
    PDF_cal12 = 0;
    for i=1:maxGen
        for j=1:ps
            if j==1 
                if pBestFit(end)<pBestFit(2), lBest=pBest(end,:);
                else, lBest=pBest(2,:);
                end
            elseif j==ps
                 if pBestFit(ps-1)<pBestFit(1),lBest=pBest(ps-1);
                 else, lBest=pBest(1,:);
                 end
            else
                 if pBestFit(j-1)<pBestFit(j+1), lBest=pBest(j-1,:);
                 else, lBest=pBest(j+1,:);
                 end
            end
            
            V(j,:)=w*V(j,:)+c1*rand(1,D).*(pBest(j,:)-pop(j,:))+c2*rand(1,D).*(lBest-pop(j,:));
            V(j,find(V(j,:)>vMax))=vMax;
            V(j,find(V(j,:)<vMin))=vMin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            fitness(j)=hybrid_func1(pop(j,:));%************************************************************************
            if fitness(j)<pBestFit(j)
                pBest(j,:)=pop(j,:);
                pBestFit(j)=fitness(j);
            end
            if fitness(j)<gBestFit
                gBest=pop(j,:);
                gBestFit=fitness(j);
            end
        end
        yy12(i+1)=gBestFit;
        if gBestFit > des_val
            PDF_cal12 = PDF_cal12 + 1;
        end     
    end
        time11(run) = toc;

    PDF_num_LPSO(run) = PDF_cal12;
    if yy12(maxGen - 5) < des_val
        PSO_12(run) = 1;
    end
    ff_end12(run,:) = yy12;
end
load('time_evoler_hybrid_L');
load('time_evoler_hybrid_D');
load('time_evoler_hybrid');

time1 = time1 + time_ini;
time2 = time2 + time_ini;
time3 = time3 + time_ini;
time4 = time4 + time_ini;
time5 = time5 + time_ini;
time6 = time6 + time_ini;
time7 = time7 + time_ini;
time8 = time8 + time_ini;
time9 = time9 + time_ini;
time10 = time10 + time_ini;
time11 = time11 + time_ini;
save time_end_all2 time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 time13 time14 time_evoler




