% clc;
clear all;
% close all;
c1 = 2;
c2 = 2;  
maxgen = 1000;       
sizepop = 20;       
D= 6;                 

Total = 1260;
num_all = xlsread('para_6.xlsx');  % load the parameter of 6 generators in the case of 6-dimensional power grid dispatch
num_test = num_all;

P_i_min_all = num_test(:,1)';
P_i_max_all = num_test(:,2)';
P_1_text = P_i_min_all(1:D-1);
P_2_text = P_i_max_all(1:D-1);
pop = [];
pop2 = [];
for iii = 1 : 50
    iii
    tran_p=0.9;
    state=1;
    neighborsize=5;
    F=0.5;
    CR=0.9;
    dF=(1-0.5)/maxgen;
    dCR=(0.4-0.9)/maxgen;
    xmax= 600*ones(1,D-1);
    xmin= 50*ones(1,D-1);
    vmax= 20;
    vmin= -20;
    x=zeros(sizepop,D-1);
    v=zeros(sizepop,D-1);
    fitness=zeros(1,sizepop);
    meandistance=zeros(1,sizepop);
    evfactor=zeros(1,sizepop);
    u=zeros(1,D-1);
    distance=zeros(1,sizepop);
    neighborindex=zeros(sizepop,neighborsize);
    status=zeros(1,sizepop);
    pbest2=zeros(1,D-1);
    gbest2=zeros(1,D-1);
    result=zeros(1,maxgen);
    for i=1:sizepop 
%         x(i,:)=(xmax-xmin).*rand(1,D)+xmin;
%         v(i,:)=rand(1,D).*(vmax-vmin)+vmin;
%         fitness(i)=rastrigin(x(i,:));%********************
        x(i,:)= (xmax-xmin).*rand(1,D-1)+xmin;
%         if i == 1
%             pop(i,:)= (x_index_all); %x_index_all(1), x_index_all(2), x_index_all(3),x_index_all(4),x_index_all(5),x_index_all(6)
%             x(i,:)= (x_index_all(1:D-1));
%         else
%             x_index_all_ini = x_index_all(1:D-1) + 4*(rand(1,D-1)*2-1);
%             pop(i,:)= ([x_index_all_ini, Total - sum(x_index_all_ini)]);           
%             x(i,:) = ([x_index_all_ini]);
%         end
        
        
        ind_1 = x(i,:) < P_1_text;
        ind_peak = find(ind_1 == 1);
        x(i,ind_peak) = P_1_text(ind_peak);
        
        ind_2 = x(i,:) > P_2_text;
        ind_peak2 = find(ind_2 == 1);
        x(i,ind_peak2) = P_2_text(ind_peak2);
        x(find(x < 0)) = 1;
        
        
        pop(i,:)= ([x(i,:), Total - sum(x(i,:))]);
        v(i,:)= rand(1,D-1).*(vmax-vmin)+vmin;
        fitness(i)= power_allocation_6D(pop(i,:));       
    end
    
    pbest=x;
    pbestfitness=fitness;
    [gbestfitness,bestindex]=min(pbestfitness);
    gbest=pbest(bestindex,:);

     for i=1:maxgen
        i;
        for j=1:sizepop
            ind_1 = x(j,:) < P_1_text;
            ind_peak = find(ind_1 == 1);
            x(j,ind_peak) = P_1_text(ind_peak);
            ind_2 = x(j,:) > P_2_text;
            ind_peak2 = find(ind_2 == 1);
            x(j,ind_peak2) = P_2_text(ind_peak2);
            x(find(x < 0)) = 1;
            pop_exi = Total - sum(x(j,:));
            pop_sum = [x(j,:), pop_exi];
            if (pop_exi >  P_i_max_all(end)) || (pop_exi < P_i_min_all(end))
                fitness(j)= 1e10;
            else
                fitness(j)= power_allocation_6D(pop_sum);
            end
            if fitness(j)<pbestfitness(j)
                pbestfitness(j)=fitness(j);
                pbest(j,:)=x(j,:);
            end
            if pbestfitness(j)<gbestfitness
                gbestfitness=pbestfitness(j);
                gbest=pbest(j,:);
            end
        end
        for j=1:sizepop
            distance=0;
            for k=1:sizepop
                distance=distance+norm(x(k,:)-x(j,:));
            end
            meandistance(j)=distance/(sizepop-1);
        end
        distancemin=min(meandistance);
        distancemax=max(meandistance);
        for j=1:sizepop
            evfactor(j)=(meandistance(j)-distancemin)/(distancemax-distancemin);
            if evfactor(j)<0.25,status(j)=1;
            elseif evfactor(j)<0.5, status(j)=2;
            elseif evfactor(j)<0.75, status(j)=3;
            else,status(j)=4;
            end
        end
        for j=1:sizepop
            for d=1:D-1
                r1=randi(sizepop);
                r2=randi(sizepop);
                u(d)=pbest(j,d)+F*(pbest(r1,d)-pbest(r2,d));
                ind_1 = u < P_1_text;
                ind_peak = find(ind_1 == 1);
                u(ind_peak) = P_1_text(ind_peak);
                ind_2 = u > P_2_text;
                ind_peak2 = find(ind_2 == 1);
                u(ind_peak2) = P_2_text(ind_peak2);
                u(find(u < 0)) = 1;
                if rand>CR, u(d)=pbest(j,d);end
            end
            
            pop_exi = Total - sum(u);
            pop_sum = [u, pop_exi];
            if (pop_exi >  P_i_max_all(end)) || (pop_exi < P_i_min_all(end))
                fitness_u= 1e10;
            else
                fitness_u= power_allocation_6D(pop_sum);
            end
%             fitness_u=rastrigin(u);%**********************************
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
        for j=1:sizepop
            for k=1:sizepop
                distance(k)=norm(pbest(j,:)-pbest(k,:));
            end
            [distance2,sortindex]=sort(distance,2);
            neighborindex(j,:)=sortindex(1:neighborsize);
        end
        for k=1:sizepop
            distance(k)=norm(gbest-pbest(k,:));
        end
        [distance2,sortindex]=sort(distance,2);
        gneighborindex=sortindex(1:neighborsize);
        for j=1:sizepop
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
        for j=1:sizepop
            if status(j)==1
                c1=2;c2=2;
                for d=1:D-1
                    pbest2(d)=pbest(j,d);
                end
                for d=1:D-1
                    gbest2(d)=gbest(d);
                end
            elseif status(j)==2
                c1=2.1;c2=1.9;
                r1=randi(sizepop);
                r2=randi(neighborsize);
                index=neighborindex(r1,r2);
                for d=1:D-1
                    pbest2(d)=pbest(index,d);
                end
                for d=1:D-1
                    gbest2(d)=gbest(d);
                end
            elseif status(j)==3
                c1=2.2;c2=1.8;
                for d=1:D-1
                    pbest2(d)=pbest(j,d);
                end
                r2=randi(neighborsize);
                index=gneighborindex(r2);
                for d=1:D-1
                    gbest2(d)=pbest(index,d);
                end
            else
                c1=1.8;c2=2.2;
                r1=randi(sizepop);
                r2=randi(neighborsize);
                index=neighborindex(r1,r2);
                for d=1:D-1
                    pbest2(d)=pbest(index,d);
                end
                r2=randi(neighborsize);
                index=gneighborindex(r2);
                for d=1:D-1
                    gbest2(d)=pbest(index,d);
                end
            end
            for d=1:D-1
                v(j,d)=w(j)*v(j,d)+c1*rand*(pbest2(d)-x(j,d))+c2*rand*(gbest2(d)-x(j,d));
                v(j,find(v(j,:)<vmin))=vmin;
                v(j,find(v(j,:)>vmax))=vmax;
                x(j,d)=x(j,d)+v(j,d);
                ind_1 = x(j,:) < P_1_text;
                ind_peak = find(ind_1 == 1);
                x(j,ind_peak) = P_1_text(ind_peak);
                ind_2 = x(j,:) > P_2_text;
                ind_peak2 = find(ind_2 == 1);
                x(j,ind_peak2) = P_2_text(ind_peak2);
                x(find(x < 0)) = 1;
                
            end
        end
        F=F+dF;
        CR=CR+dCR;
     end
%     g_Best_end = [gBest, (Total - sum(gBest))];
    Y(iii,:) = result;
%     min(result)
%     hold on
%     plot(result,'r')
%     xlabel('Genration');
%     ylabel('Fitness')
end
save('DNSPSO_6D','Y')

