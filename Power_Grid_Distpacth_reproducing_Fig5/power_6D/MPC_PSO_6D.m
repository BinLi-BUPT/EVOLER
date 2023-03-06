% clc;
clear all;
% close all;
maxgen = 1000;          
sizepop = 20;       
D = 6;                 
Total = 1260;
num_all = xlsread('para_6.xlsx'); % load the parameters of 6 power generators
num_test = num_all;

P_i_min_all = num_test(:,1)';
P_i_max_all = num_test(:,2)';
P_1_text = P_i_min_all(1:D-1);
P_2_text = P_i_max_all(1:D-1);
pop = [];
pop2 = [];

for iii = 1 : 50
    iii
    t=0.5;
    MAXT=5;
    F=0.5;
    N1=t*sizepop;%
    N2=sizepop-N1;
    c=1.49445;
    w=0.729;
    tournamentsize=5;
    popmax= P_2_text;
    popmin= P_1_text;
    vmax=0.2*popmax;
    vmin=0.2*popmin;
    pop=zeros(sizepop,D-1);
    v=zeros(sizepop,D-1);
    fitness=zeros(1,sizepop);
    pbestfitness=zeros(1,sizepop);
    pbest=zeros(sizepop,D-1);
    gbest=zeros(1,D-1);
    mbest=zeros(1,D-1);
    result=zeros(1,maxgen);
    
    for i=1:sizepop
        while(1)
            for d=1:D-1
                pop(i,d)=rand*(popmax(d)-popmin(d))+popmin(d);
                v(i,d)=rand*(vmax(d)-vmin(d))+vmin(d);
                %pop(i,d)=round(pop(i,d));
            end
            ind_1 = pop(i,:) < P_1_text;
            ind_peak = find(ind_1 == 1);
            pop(i,ind_peak) = P_1_text(ind_peak);

            ind_2 = pop(i,:) > P_2_text;
            ind_peak2 = find(ind_2 == 1);
            pop(i,ind_peak2) = P_2_text(ind_peak2);
            pop(find(pop < 0)) = 1;
            if Total - sum(pop(i,:))>=P_i_min_all(end) && Total - sum(pop(i,:))<=P_i_max_all(end), break; end
        end
        pop3(i,:)= ([pop(i,:), Total - sum(pop(i,:))]);
        fitness(i)= power_allocation_6D(pop3(i,:));       
    end
    me1=sum(fitness);
    pbestfitness=fitness;
    pbest=pop;
    [sortfitness,sortindex]=sort(pbestfitness,2);
    gbestfitness=sortfitness(1);
    gbest=pbest(sortindex(1),:);
    T=0;
    
    
    %%
    for i=1:maxgen
        i;
        is_des=0;%
        [sortfitness,sortindex]=sort(pbestfitness,2);
        EP=sortindex(1:N1);
        GP=sortindex(N1+1:end);
        segmentnum=ceil((i/maxgen)*(D-1));
        ds=ceil((D-1)/segmentnum);
        y=floor((D-1)/segmentnum);
        sumfitness=0;
        for j=1:N1
            sumfitness=sumfitness+pbestfitness(GP(j));
        end
        phi=-1*me1*sumfitness/N1;
        phi=1/((1+exp(phi))^i);%
        for j=1:N2
            index=GP(j);
            for dd=1:ds:(D-1)
                sum1=0;
                sum2=zeros(1,D-1);
                pl=zeros(1,D-1);
                for ii=1:tournamentsize
                    tournamentindex=EP(randi(N1));
                    sum1=sum1+pbest(tournamentindex);
                    sum2=sum2+pbestfitness(tournamentindex)*pbest(tournamentindex,:);
                end
                pl=sum2/(sum1*tournamentsize);
                r1=rand;r2=rand;
                if dd+ds-1>D-1
                    for d=(D-1-y+1):D-1
                        v(index,d)=w*v(index,d)+c*r1*(pl(d)-pop(index,d))+phi*r2*(gbest(d)-pop(index,d));
                        if v(index,d)<vmin(d), v(index,d)=vmin(d); end
                        if v(index,d)>vmax(d), v(index,d)=vmax(d); end
                        pop(index,d)=pop(index,d)+v(index,d);
                        if pop(index,d)<popmin(d), pop(index,d)=popmin(d); end
                        if pop(index,d)>popmax(d), pop(index,d)=popmax(d); end
                    end
                else
                    for d=dd:(dd+ds-1)
                        v(index,d)=w*v(index,d)+c*r1*(pl(d)-pop(index,d))+phi*r2*(gbest(d)-pop(index,d));
                        if v(index,d)<vmin(d), v(index,d)=vmin(d); end
                        if v(index,d)>vmax(d), v(index,d)=vmax(d); end
                        pop(index,d)=pop(index,d)+v(index,d);
                        if pop(index,d)<popmin(d), pop(index,d)=popmin(d); end
                        if pop(index,d)>popmax(d), pop(index,d)=popmax(d); end
                    end
                end
            end
            
            pop3_temp= ([pop(index,:), Total - sum(pop(index,:))]);
            fitness(index)= power_allocation_6D(pop3_temp);
            
            
            %  fitness(index)=griewank(pop(index,:));%***************************************************************
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
            phi=1/((1+exp(phi))^i);%
            sumpop=zeros(1,D-1);
            for j=1:N1
                index=EP(j);
                sumpop=sumpop+pbest(index,:);
            end
            for j=1:N1
                index=EP(j);
                sumD=sum(pbest(index,:));
                for d=1:D-1
                    lambda=1/(exp(sumD/(D-1)-pbest(index,d))+1);
                    r=rand;
                    mbest(d)=lambda*r/(D-1)*sumD+(1-lambda)*(1-r)*sumpop(d);
                    v(index,d)=w*v(index,d)+c*rand*(mbest(d)-pop(index,d));
                    if v(index,d)<vmin(d), v(index,d)=vmin(d); end
                    if v(index,d)>vmax(d), v(index,d)=vmax(d); end
                    pop(index,d)=phi*pop(index,d)+v(index,d);
                    if pop(index,d)<popmin(d), pop(index,d)=popmin(d); end
                    if pop(index,d)>popmax(d), pop(index,d)=popmax(d); end
                end
                pop3_temp= ([pop(index,:), Total - sum(pop(index,:))]);
                fitness(index)= power_allocation_6D(pop3_temp);       %
                
                %fitness(index)=griewank(pop(index,:));%***************************************************************
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
                index1=randi(sizepop);
                index2=randi(sizepop);
                while(index1==index || index2==index || index1==index2)
                    index1=randi(sizepop);
                    index2=randi(sizepop);
                end
                for d=1:(D-1)
                    pop(index,d)=gbest(d)+F*(pbest(index1,d)-pbest(index2,d));
                    if pop(index,d)<popmin(d), pop(index,d)=popmin(d); end
                    if pop(index,d)>popmax(d), pop(index,d)=popmax(d); end
                end
                pop3_temp= ([pop(index,:), Total - sum(pop(index,:))]);
                fitness(index)= power_allocation_6D(pop3_temp);       %
                % fitness(index)=griewank(pop(index,:));%********************************************************************
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
        T;
    end
    
    Y(iii,:) = result;
    %     hold on
    %     plot(result,'r')
    %     xlabel('Genration');
    %     ylabel('Fitness')
end
% Y_2 = mean(Y,1);
% figure
% plot(Y_2,'r')
% xlabel('Genration');
% ylabel('Fitness')

save('MPCPSO_6D','Y')

% de1 = load('GC_PSO_data')