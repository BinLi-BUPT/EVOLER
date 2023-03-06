% Copyright 2022, All Rights Reserved
% Code by Bin Li, Ziping Wei

%% Several Comparing methods
%Standard PSO、CLPSO、CLPSO-H、CLPSO-S、ML-PSO、DPSO、NPSO、GCPSO、MPCPSO、DNSPSO、Local PSO
clc;
clear all;
close all
%% 
N_add=4;
des_val= 1e-19;   %If the output result is less than the threshold 'des_val', the global optimum is found  
popMax= 5;       %Variable Range [popmin: popmax]

D=2;
ps=50;
maxGen=500;
num_all = 500;
popMin= - popMax; %Variable Range [popmin: popmax]
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


for run=1:num_all
    run
    M= (-1+2*rand(1,D));
    for i=1:ps
        initPop(i,:)=popMin+(popMax-popMin)*rand(1,D);
        initV(i,:)=vMin+(vMax-vMin)*rand(1,D);
        initFitness(i)=shifted_rastrigin(initPop(i,:),M);%*****************************************************************
    end
    initpBest=initPop;
    initpBestFit=initFitness;
    [initgBestFit,index]=min(initpBestFit);
    initgBest=initpBest(index);
    %% Standard PSO 1
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
            fitness(j)=shifted_rastrigin(pop(j,:),M);%************************************************************************
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
    PDF_num_SPSO(run) = PDF_cal;
    if yy1(maxGen - 5) < des_val
        PSO_1(run) = 1;
    end
    ff_end1(run,:) = yy1;
    %% DPSO 2
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
            fitness(j)=shifted_rastrigin(pop(j,:),M);%**********************************************************************
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
    PDF_num_DPSO(run) = PDF_cal2;
    if yy2(maxGen - 5) < des_val
        PSO_2(run) = 1;
    end
    ff_end2(run,:) = yy2;
    %% CPSO-S 3
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
                 fit1=shifted_rastrigin(temp1,M);
                 fit2=shifted_rastrigin(pBest1,M);%***************************************************************
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
     %% CPSO-H 4
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
                fit1=shifted_rastrigin(temp1,M);fit2=shifted_rastrigin(pBest1,M);%****************************************************
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
            s_fitness(i)=shifted_rastrigin(s_pop(i,:),M);%*************************************************************
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
    PDF_num_CPSOH(run) = PDF_cal4;
    if yy4(maxGen - 5) < des_val
        PSO_4(run) = 1;        
    end
    ff_end4(run,:) = yy4;
    %% CLPSO 5
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
            fitness(i)=shifted_rastrigin(pop(i,:),M);%****************************************************************
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

    if yy5(maxGen - 5) < des_val
        PSO_5(run) = 1;
    end
    ff_end5(run,:) = yy5;
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
                Fitness(j,k)=shifted_rastrigin(reshape(Pop(j,k,:),1,D),M);%******************************************************************************

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
            fit=shifted_rastrigin(reshape(Pop(groupNum,k,:),1,D),M);%*************************************************************************************
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
    PDF_num_MLPSO(run) = PDF_cal6;
    if yy6(maxGen - 5) < des_val
        PSO_6(run) = 1;
    end
    ff_end6(run,:) = yy6;
    %% NPSO 7
    c1=2;c2=2;
    pop=initPop;  
    V =initV;
    fitness=initFitness;
    fitnesspbest=fitness;
    gbest=initgBest;
    fitnessgbest=initgBestFit;
    [fitnessgworst, bestindex2]=max(initFitness);
    g_worst=initPop(bestindex2,:);
    p_worst=initPop;
    fitnesspworst=fitness;
    
    PDF_cal7 = 0;
    yy7(1)=fitnessgbest;
    for i=1:maxGen
        for j=1:ps
            V(j,:)=V(j,:)+c1*rand(1,D).*(pop(j,:)-p_worst(j,:))+c2*rand(1,D).*(pop(j,:)-g_worst);
            V(j,find(V(j,:)>vMax))=vMax;
            V(j,find(V(j,:)<vMin))=vMin;
            pop(j,:) = pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
            
            fitness(j)=shifted_rastrigin(pop(j,:),M);
            if fitness(j) < fitnesspbest(j)
                fitnesspbest(j)=fitness(j);
                pbest(j,:)=pop(j,:);
%             else
%                 p_worst(j,:) = pop(j,:);      
            end
            if fitness(j)>fitnesspworst(j)
                fitnesspworst(j)=fitness(j);
                p_worst(j,:)=pop(j,:);
            end
            if fitness(j) < fitnessgbest
                fitnessgbest=fitness(j);
                gbest=pop(j,:);
%             else
%                 g_worst = pop(j,:);
            end  
            if fitness(j)>fitnessgworst
                fitnessgworst=fitness(j);
                g_worst=pop(j,:);
            end
        end
        yy7(i+1)=fitnessgbest;
        if fitnessgbest > des_val
            PDF_cal7 = PDF_cal7 + 1;
        end
    end  %%
    PDF_num_NPSO(run) = PDF_cal7;
    if yy7(maxGen - 5) < des_val
        PSO_7(run) = 1;
    end
    ff_end7(run,:) = yy7;
    %% GCPSO 8
    
    c1=1.49; c2=1.49; w=0.72;
    ft = 1;                                 %range control
    success = 15;                           %numbers of success
    failure = 5;                            %numbers of failure
    suNumber = 0;                           %count the number of success
    faNumber = 0;                           %count the number of failure
    speedKeep = zeros(1,D);        %keep the speed of best particle
    fitness=initFitness;
    pop=initPop;
    gbest=initgBest;
    fitnessgbest=initgBestFit;
    speed =initV;
    pbest=pop;
    fitnesspbest=fitness;
    nMin = index;
    %initial data storage
    xBest = initPop;                                 %each particle's best position
    yBest = zeros(1,D);                              %best position for each particles in history
    xFit =  fitnesspbest;                            %each particle's best fittness
    yFit =  fitnessgbest;                            %best fittness in history
    yy8(1) = initgBestFit;
    PDF_cal8 = 0;
    
    for iIteration = 1:maxGen
        for iParticles = 1:ps
            %xFitness = shifted_rastrigin();
            xFitness=shifted_rastrigin(pop(iParticles,:),M);    
            %current particle fittness
            if xFitness < xFit(iParticles)            % renew particle's best position
                xFit(iParticles) = xFitness;
                xBest(iParticles,:) = pop(iParticles,:);
            end
        end
        
        if min(xFit)  < yFit                          % renew best fittness
            [yFit,nMin] = min(xFit);
            yBest = xBest(nMin,:);
            %%first insert
            suNumber = suNumber + 1;
            faNumber = 0;
        else
%             [yFit,nMin] = min(xFit);
%             yBest = xBest(nMin,:);
            faNumber = faNumber + 1;
            suNumber = 0;
        end
        %%iteration end,renew speed and position
        %second insert
        speed = (speed * w + c1 * rand * (xBest - pop) + c2 * rand * (repmat(yBest,ps,1) - pop));%renew speed
        %boundary speed process
        speed(speed < vMin ) = vMin;
        speed(speed > vMax) = vMax;
        speedKeep = speed(nMin,:);
        pop = pop + speed;          %renew position
        %third insert
        if suNumber > success
            ft = 2 * ft;
            failure = failure + 1;  %renew failure 
        elseif faNumber > failure
            ft = 0.5 * ft;
            success = success + 1;  %renew success 
        end
        pop(nMin,:) = yBest(1,:) + w * speedKeep + ft .* (1 - 2 * rand(1,D));
        for j=1:ps
            pop(j,find(pop(j,:)>popMax))=popMax;
            pop(j,find(pop(j,:)<popMin))=popMin;
        end
        yy8(iIteration+1)=yFit;
%         if fitnessgbest > des_val
        if yFit>des_val
            PDF_cal8 = PDF_cal8 + 1;
        end
    end  %
    PDF_num_GCPSO(run) = PDF_cal8;
    if yy8(maxGen - 5) < des_val
        PSO_8(run) = 1;
    end
    ff_end8(run,:) = yy8;
    %% MPCPSO 9
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
    yy9(1)=gbestfitness;
    PDF_cal9 = 0;
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
            fitness(index)=shifted_rastrigin(pop(index,:),M);%*****************************************************************
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
                 fitness(index)=shifted_rastrigin(pop(index,:),M);%*****************************************************************
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
                 fitness(index)=shifted_rastrigin(pop(index,:),M);%*****************************************************************
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
        yy9(i+1)=gbestfitness;
        if gbestfitness > des_val
            PDF_cal9 = PDF_cal9 + 1;
        end
    end
    PDF_num_MPCPSO(run) = PDF_cal9;
     if yy9(maxGen - 5) < des_val
        PSO_9(run) = 1;
    end
    ff_end9(run,:) = yy9;
    %% DNSPSO 10
    PDF_cal10 = 0;
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
    yy10(1)=gbestfitness;
    for i=1:maxGen
        i;
        for j=1:ps
            fitness(j)=shifted_rastrigin(pop(j,:),M);%***********************
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
            fitness_u=shifted_rastrigin(u,M);%**********************************
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
        yy10(i+1)=gbestfitness;
        if gbestfitness > des_val
            PDF_cal10 = PDF_cal10 + 1;
        end
    end
    PDF_num_DNSPSO(run) = PDF_cal10;
    if yy10(maxGen - 5) < des_val
        PSO_10(run) = 1;
    end
    ff_end10(run,:) = yy10;
    t11=clock;
    t11=t11(4:6);
    %% Local PSO 11
    c1=1.4962; c2=1.4962; w=0.7968;
    pop=initPop;
    V=initV;
    fitness=initFitness;
    pBest=pop;
    pBestFit=fitness;
    [gBestFit,index]=min(pBestFit);
    gBest=pBest(index);
    yy11(1)=gBestFit;
    PDF_cal11 = 0;
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
            fitness(j)=shifted_rastrigin(pop(j,:),M);%************************************************************************
            if fitness(j)<pBestFit(j)
                pBest(j,:)=pop(j,:);
                pBestFit(j)=fitness(j);
            end
            if fitness(j)<gBestFit
                gBest=pop(j,:);
                gBestFit=fitness(j);
            end
        end
        yy11(i+1)=gBestFit;
        if gBestFit > des_val
            PDF_cal11 = PDF_cal11 + 1;
        end     
    end
    PDF_num_LPSO(run) = PDF_cal11;
    if yy11(maxGen - 5) < des_val
        PSO_11(run) = 1;
    end
    ff_end11(run,:) = yy11;
    
end
shifted_rastrigin_2D_fitness={};


shifted_rastrigin_2D_fitness{1}=ff_end1;%standard pso
shifted_rastrigin_2D_fitness{2}=ff_end2;%dpso
shifted_rastrigin_2D_fitness{3}=ff_end3;%cpso-s
shifted_rastrigin_2D_fitness{4}=ff_end4;%cpso-h
shifted_rastrigin_2D_fitness{5}=ff_end5;%clpso
shifted_rastrigin_2D_fitness{6}=ff_end6;%mcpso
shifted_rastrigin_2D_fitness{7}=ff_end7;%npso
shifted_rastrigin_2D_fitness{8}=ff_end8;%gcpso
shifted_rastrigin_2D_fitness{9}=ff_end9;%mpcpso
shifted_rastrigin_2D_fitness{10}=ff_end10;%dnspso
shifted_rastrigin_2D_fitness{11}=ff_end11;%local pso

da_1_Pro = load('Pro_data');
ff_end12=da_1_Pro.Fit_pro(1:num_all,:);
shifted_rastrigin_2D_fitness{12}=ff_end12;%pso evoler
da_l_Pro2=load('Pro_Local_data');
ff_end13=da_l_Pro2.Fit_pro(1:num_all,:);
shifted_rastrigin_2D_fitness{13}=ff_end13;%local pso evoler
da_l_Pro3=load('Pro_downhill_data');
ff_end14=da_l_Pro3.Fit_pro(1:num_all,:);
shifted_rastrigin_2D_fitness{14}=ff_end14;
save('shifted_rastrigin_2D_fitness','shifted_rastrigin_2D_fitness')
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'des_val'},'2D_shif_rast_fitness','A1');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',des_val,'2D_shif_rast_fitness','B1');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'method'},'2D_shif_rast_fitness','A3');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'fitness'},'2D_shif_rast_fitness','B3');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'standard PSO'},'2D_shif_rast_fitness','A5');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{1},'2D_shif_rast_fitness','B5');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'DPSO'},'2D_shif_rast_fitness','A510');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{2},'2D_shif_rast_fitness','B510');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'CPSO-S'},'2D_shif_rast_fitness','A1015');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{3},'2D_shif_rast_fitness','B1015');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'CPSO-H'},'2D_shif_rast_fitness','A1520');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{4},'2D_shif_rast_fitness','B1520');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'CLPSO'},'2D_shif_rast_fitness','A2025');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{5},'2D_shif_rast_fitness','B2025');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'MCPSO'},'2D_shif_rast_fitness','A2530');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{6},'2D_shif_rast_fitness','B2530');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'NPSO'},'2D_shif_rast_fitness','A3035');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{7},'2D_shif_rast_fitness','B3035');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'GCPSO'},'2D_shif_rast_fitness','A3540');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{8},'2D_shif_rast_fitness','B3540');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'MPCPSO'},'2D_shif_rast_fitness','A4045');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{9},'2D_shif_rast_fitness','B4045');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'DNSPSO'},'2D_shif_rast_fitness','A4550');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{10},'2D_shif_rast_fitness','B4550');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'Local PSO'},'2D_shif_rast_fitness','A5055');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{11},'2D_shif_rast_fitness','B5055');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'EVOLER'},'2D_shif_rast_fitness','A5560');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{12},'2D_shif_rast_fitness','B5560');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'EVOLER_Local'},'2D_shif_rast_fitness','A6065');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{13},'2D_shif_rast_fitness','B6065');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',{'EVOLER_Downhill'},'2D_shif_rast_fitness','A6570');
s=xlswrite('shifted_rastrigin_2D_fitness.xlsx',shifted_rastrigin_2D_fitness{14},'2D_shif_rast_fitness','B6570');



pro_12 = mean(ff_end12);
PSO_12 = [];
for ii = 1 :size(ff_end12,1)
    PDF_cal12 = 0;
    dp_now = ff_end12(ii,:);
    for jj = 1 : length(dp_now)
        if dp_now(jj) > des_val
            PDF_cal12 = PDF_cal12 + 1;
        end
    end
    PDF_num_Pro12(ii) = PDF_cal12;
    if ff_end12(ii,end-3) < des_val
        PSO_12(ii) = 1;
    end
end


pro_13 = mean(ff_end13);
PSO_13 = [];
for ii = 1 :size(ff_end13,1)
    PDF_cal13 = 0;
    dp_now = ff_end13(ii,:);
    for jj = 1 : length(dp_now)
        if dp_now(jj) > des_val
            PDF_cal13 = PDF_cal13 + 1;
        end
    end
    PDF_num_Pro13(ii) = PDF_cal13;
    if ff_end13(ii,end-3) < des_val
        PSO_13(ii) = 1;
    end
end
pro_14 = mean(ff_end14);
PSO_14 = [];
for ii = 1 :size(ff_end14,1)
    PDF_cal14 = 0;
    dp_now = ff_end14(ii,:);
    for jj = 1 : length(dp_now)
        if dp_now(jj) > des_val
            PDF_cal14 = PDF_cal14 + 1;
        end
    end
    PDF_num_Pro14(ii) = PDF_cal14;
    if ff_end14(ii,end-3) < des_val
        PSO_14(ii) = 1;
    end
end
% X1 = PDF_num_Pro;

figure()
pop=1:maxGen+1;
xlim([0,900]);
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end1,1),'d-','Linewidth',1,'Color',[1.00,0.41,0.16],'MarkerIndices',marker_idx);%standard pso
hold on;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end2,1),'v-','Linewidth',1,'Color',[0.08,0.71,0.68],'MarkerIndices',marker_idx);%dpso
hold on;
marker_idx=1:30:maxGen*ps*2;
semilogy(pop*ps*2,mean(ff_end3,1),'h-','Linewidth',1,'Color',[0.06,0.97,0.60],'MarkerIndices',marker_idx);%cpso-s
hold on;
marker_idx=1:30:maxGen*ps*1.5;
semilogy(pop*ps*1.5,mean(ff_end4,1),'x-','Linewidth',1,'Color',[0.57,0.92,0.57],'MarkerIndices',marker_idx);%cpso-h
hold on;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end5,1),'+-','Linewidth',1,'Color',[0.50,0.50,0.50],'MarkerIndices',marker_idx);%clpso
hold on;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end6,1),'*-','Linewidth',1,'Color',[0.93,0.69,0.13],'MarkerIndices',marker_idx);%mcpso
hold on;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end7,1),'.-','Linewidth',1.5,'Color',[0.96,0.46,0.77],'MarkerIndices',marker_idx)%npso
hold on
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps, mean(ff_end8,1),'v-','Linewidth',1,'Color',[0.21,0.21,0.85],'MarkerIndices',marker_idx);%gcpso
hold on;
marker_idx=1:30:maxGen*ps*2;
semilogy(pop*ps*2, mean(ff_end9,1),'o-','Linewidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',marker_idx);%mpcdns
hold on
marker_idx=1:30:maxGen*ps*2;
semilogy(pop*ps*2, mean(ff_end10,1),'--','Linewidth',2,'Color',[0.72,0.27,1.00],'MarkerIndices',marker_idx);%dnspso
hold on;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end11,1),'-p','Linewidth',1,'Color',[0.94,0.50,0.94],'MarkerIndices',marker_idx);%lpso
hold on;

pop=1:maxGen;
marker_idx=1:30:maxGen*ps;
semilogy(pop*ps,mean(ff_end12,1),'-s','Linewidth',1,'Color','r','MarkerIndices',marker_idx);%evoler
hold on;
semilogy(pop*ps,mean(ff_end13,1),'-.','Linewidth',1,'Color',[1.00,0.28,0.75],'MarkerIndices',marker_idx);%local pso evoler
hold on;

Xdata_de=1:maxGen+1;
Xdata_de(1:N_add)=Xdata_de(1:N_add)*ps;
Xdata_de(N_add+1:501)=Xdata_de(N_add)+(Xdata_de(N_add+1:end)-N_add)*2;
semilogy(Xdata_de,mean(ff_end14,1),'Color',[0.90,0.63,0.75],'Linewidth',1,'MarkerIndices',marker_idx);%downhill evoler
hold on;

xlabel('Number of function evaluations');
ylabel('Fitness');
legend('Standard PSO','DPSO','CPSO-S','CPSO-H','CLPSO','MCPSO','NPSO','GCPSO','MPCPSO','DNSPSO','Local PSO','EVOLER','EVOLER_l','EVOLER_d','location','SouthWest');
xlim([0,maxGen*ps]);


[~,id_1] = find(PSO_1 == 1);
[~,id_2] = find(PSO_2 == 1);
[~,id_3] = find(PSO_3 == 1);
[~,id_4] = find(PSO_4 == 1);
[~,id_5] = find(PSO_5 == 1);
[~,id_6] = find(PSO_6 == 1);
[~,id_7] = find(PSO_7 == 1);
[~,id_8] = find(PSO_8 == 1);
[~,id_9] = find(PSO_9 == 1);
[~,id_10] = find(PSO_10 == 1);
[~,id_11] = find(PSO_11 == 1);
[~,id_12] = find(PSO_12 == 1);
[~,id_13] = find(PSO_13 == 1);
[~,id_14] = find(PSO_14 == 1);

PDF_num_SPSO_1 = PDF_num_SPSO(id_1);
[y_pso,x_pso]=hist(PDF_num_SPSO_1,20);
if isempty(PDF_num_SPSO_1)
    PDF_num_SPSO_1 = maxGen;
    p1 = 0;
else
    p1 = length(PDF_num_SPSO_1)/num_all;
end

PDF_num_DPSO_2 = PDF_num_DPSO(id_2);
[y_dpso,x_dpso]=hist(PDF_num_DPSO_2,20);
if isempty(PDF_num_DPSO_2)
    PDF_num_DPSO_2 = maxGen;
    p2 = 0;
else
    p2 = length(PDF_num_DPSO_2)/num_all;
end


PDF_num_CPSO_3 = PDF_num_CPSO(id_3);
[y_cpsos,x_cpsos]=hist(PDF_num_CPSO_3,20);
if isempty(PDF_num_CPSO_3)
    PDF_num_CPSO_3 = maxGen;
    p3 = 0;
else
    p3 = length(PDF_num_CPSO_3)/num_all;
end


PDF_num_CPSOH_4 = PDF_num_CPSOH(id_4);
[y_cpsoh,x_cpsoh]=hist(PDF_num_CPSOH_4,20);
if isempty(PDF_num_CPSOH_4)
    PDF_num_CPSOH_4 = maxGen;
    p4 = 0;
else
    p4 = length(PDF_num_CPSOH_4)/num_all;
end


PDF_num_CLPSO_5 = PDF_num_CLPSO(id_5);
[y_clpso,x_clpso]=hist(PDF_num_CLPSO_5,20);
if isempty(PDF_num_CLPSO_5)
    PDF_num_CLPSO_5 = maxGen;
    p5 = 0;
else
    p5 = length(PDF_num_CLPSO_5)/num_all;
end

PDF_num_MLPSO_6 = PDF_num_MLPSO(id_6);
[y_mlpso,x_mlpso]=hist(PDF_num_MLPSO_6,20);
if isempty(PDF_num_MLPSO_6)
    PDF_num_MLPSO_6 = maxGen;
    p6 = 0;
else
    p6 = length(PDF_num_MLPSO_6)/num_all;
end

PDF_num_NPSO_7 = PDF_num_NPSO(id_7);
[y_nppso,x_nppso]=hist(PDF_num_NPSO_7,20);
if isempty(PDF_num_NPSO_7)
    PDF_num_NPSO_7 = maxGen;
    p7 = 0;
else
    p7 = length(PDF_num_NPSO_7)/num_all;
end

PDF_num_GCPSO_8 = PDF_num_GCPSO(id_8);
[y_gcpso,x_gcpso]=hist(PDF_num_GCPSO_8,20);
if isempty(PDF_num_GCPSO_8)
    PDF_num_GCPSO_8 = maxGen;
    p8 = 0;
else
    p8 = length(PDF_num_GCPSO_8)/num_all;
end

PDF_num_MPCPSO_9 = PDF_num_MPCPSO(id_9);
[y_pso,x_pso]=hist(PDF_num_MPCPSO_9,20);
if isempty(PDF_num_MPCPSO_9)
    PDF_num_MPCPSO_9 = maxGen;
    p9 = 0;
else
    p9 = length(PDF_num_MPCPSO_9)/num_all;
end

PDF_num_DNSPSO_10 = PDF_num_DNSPSO(id_10);
[y_pso,x_pso]=hist(PDF_num_DNSPSO_10,20);
if isempty(PDF_num_DNSPSO_10)
    PDF_num_DNSPSO_10 = maxGen;
    p10 = 0;
else
    p10 = length(PDF_num_DNSPSO_10)/num_all;
end

PDF_num_LPSO_11 = PDF_num_LPSO(id_11);
[y_pso,x_pso]=hist(PDF_num_LPSO_11,20);
if isempty(PDF_num_LPSO_11)
    PDF_num_LPSO_11 = maxGen;
    p11 = 0;
else
    p11 = length(PDF_num_LPSO_11)/num_all;
end


PDF_num_EVOLER_12 = PDF_num_Pro12(id_12);
[y_pso,x_pso]=hist(PDF_num_EVOLER_12,20);
if isempty(PDF_num_EVOLER_12)
    PDF_num_EVOLER_12 = maxGen;
    p12 = 0;
else
    p12 = length(PDF_num_EVOLER_12)/num_all;
end

PDF_num_EVOLER_13 = PDF_num_Pro13(id_13);
[y_pso,x_pso]=hist(PDF_num_EVOLER_13,20);
if isempty(PDF_num_EVOLER_13)
    PDF_num_EVOLER_13 = maxGen;
    p13 = 0;
else
    p13 = length(PDF_num_EVOLER_13)/num_all;
end

PDF_num_EVOLER_14 = PDF_num_Pro14(id_14);
[y_pso,x_pso]=hist(PDF_num_EVOLER_14,20);
if isempty(PDF_num_EVOLER_14)
    PDF_num_EVOLER_14 = maxGen;
    p14 = 0;
else
    p14= length(PDF_num_EVOLER_14)/num_all;
end



mean_various = [mean(PDF_num_SPSO_1),mean(PDF_num_DPSO_2),mean(PDF_num_CPSO_3),mean(PDF_num_CPSOH_4),mean(PDF_num_CLPSO_5),mean(PDF_num_MLPSO_6),mean(PDF_num_NPSO_7),mean(PDF_num_GCPSO_8),mean(PDF_num_MPCPSO_9),mean(PDF_num_DNSPSO_10),mean(PDF_num_LPSO_11),mean(PDF_num_EVOLER_12),mean(PDF_num_EVOLER_13),mean(PDF_num_EVOLER_14)];
Probablity = [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14];
figure(2)
h1 = plot(mean_various(1),Probablity(1),'x','MarkerSize',12,'lineWidth',1.5,'color',[1,0.411764705882353,0.160784313725490]);%standard pso
set(h1,'MarkerFaceColor',get(h1,'color'));
hold on
h2 = plot(mean_various(2),Probablity(2),'+','MarkerSize',12,'lineWidth',2,'color',[1,0,0]);%dpso
set(h2,'MarkerFaceColor',get(h2,'color'));
hold on
h3 = plot(mean_various(3),Probablity(3),'square','MarkerSize',12,'lineWidth',2,'color',[0.521568627450980,1,0.168627450980392]);%cpso-s
set(h3,'MarkerFaceColor',get(h3,'color'));
hold on
h4 = plot(mean_various(4),Probablity(4),'diamond','MarkerSize',12,'lineWidth',2,'color',[0.250980392156863,0.780392156862745,0.250980392156863]);%dpso-h
set(h4,'MarkerFaceColor',get(h4,'color'));
hold on
h5 = plot(mean_various(5),Probablity(5),'square','MarkerSize',8,'lineWidth',1.5,'color',[0.149019607843137,0.149019607843137,0.149019607843137]);%clpso
set(h5,'MarkerFaceColor',get(h5,'color'));
hold on
h6 = plot(mean_various(6),Probablity(6),'*','MarkerSize',8,'lineWidth',1.5,'color',[0.929411764705882,0.694117647058824,0.125490196078431]);%mcpso
set(h6,'MarkerFaceColor',get(h6,'color'));
hold on
h7 = plot(mean_various(7),Probablity(7),'^','MarkerSize',8,'lineWidth',1.5,'color',[0.929411764705882,0.450980392156863,0.749019607843137]);%npso
set(h7,'MarkerFaceColor',get(h7,'color'));
hold on
h8 = plot(mean_various(8),Probablity(8),'hexagram','MarkerSize',12,'lineWidth',1.5,'color',[0,0,1]);%gcpso
set(h8,'MarkerFaceColor',get(h8,'color'));
hold on;
h9 = plot(mean_various(9),Probablity(9),'>','MarkerSize',12,'lineWidth',1.5,'color',[0.6,0.411764705882353,0.160784313725490]);%mpcpso
set(h9,'MarkerFaceColor',get(h9,'color'));
hold on
h10 = plot(mean_various(10),Probablity(10),'<','MarkerSize',12,'lineWidth',1.5,'color',[1,0.411764705882353,0.160784313725490]);%dnspso
set(h10,'MarkerFaceColor',get(h10,'color'));
hold on
h11 = plot(mean_various(11),Probablity(11),'*','MarkerSize',12,'lineWidth',1.5,'color',[0,1,0]);%lpso
set(h11,'MarkerFaceColor',get(h11,'color'));
hold on;
h12 = plot(mean_various(12),Probablity(12),'o','MarkerSize',10,'lineWidth',1.5,...
                'color',[0.0784313725490200,0.709803921568627,0.678431372549020]);%pso evoler
set(h12,'MarkerFaceColor',get(h12,'color'));
hold on;
h13 = plot(mean_various(13),Probablity(13),'x','MarkerSize',12,'lineWidth',2,'color',[0.15,0.15,0.15]);%local pso evoler
set(h13,'MarkerFaceColor',get(h13,'color'));
hold on;
h14 = plot(mean_various(14),Probablity(14),'>','MarkerSize',12,'lineWidth',1.5,'color',[0.70,0.50,0.70]);%downhill evoler 
set(h14,'MarkerFaceColor',get(h14,'color'));
hold on;

legend('Standard PSO','DPSO','CPSO-S','CPSO-H','CLPSO','MCPSO','NPSO','GCPSO','MPCPSO','DNSPSO','Local PSO','EVOLER','EVOLER_l','EVOLER_d','location','SouthWest');


