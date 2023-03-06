% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clear all;
close all;

%-----------initialization parameter setting-------------------------
maxg = 200;                  % maximum number of generations
c1 = linspace(2,2.1,maxg);   % Cognition learning factor c1
c2 = linspace(1.8,1.9,maxg); % Social learning factor c2
sizepop = 20;                % population size
D = 6;                       % dimension
Vmax = 10;
Vmin = -10;
Total = 1260;

%(i): load the estimated global optimum in the representation space. Note that, the low-rank represention of original problem space
% has been done by the source code file 'Low_Rank_Represention'.
load('6D_result');     % run the file 'Low_Rank_Represention' to generate this data '6D_result'
x_index_all = P_find;  % the center of the attention subspace, i.e., the global optimum in the representation space
pop = [];
pop2 = [];
N_all=  50;            % the iteration times for performance evaluate
N_add = floor((144*2*2+253*2*2+128*2*2-2*8)/sizepop); % the overall required samples to reconstruct an attention subspace
Y = 1.58e4*ones(N_all,N_add+maxg);

%(ii): Evolutionary PSO method, which is used to exploit the identified attention subspace and finally gain the global optimum.
% Initialization of population, around the identified attention subspace
for iii = 1 : N_all
    iii
    for i=1:sizepop 
        if i == 1
            pop(i,:)= (x_index_all);  
            pop2(i,:)= (x_index_all(1:D-1));
        else
            x_index_all_ini = x_index_all(1:D-1) + 4*(rand(1,D-1)*2-1);
            pop(i,:)= ([x_index_all_ini, Total - sum(x_index_all_ini)]);         
            pop2(i,:) = ([x_index_all_ini]);
        end
        V(i,:)=(rand(1,D-1)*2-1);               
        fitness(i)= power_allocation_6D(pop(i,:));      
    end
    pBest= pop2;   
    fitnesspbest=fitness;      
    
    [fitnessgbest bestindex]=min(fitness);
    gBest=(pop2(bestindex,:));
    result = zeros(1,maxg);
    P_i_min_all = [100;50;80;50;50;50]';      %lower power  bound of each generator
    P_i_max_all = [500;200;300;150;200;120]'; %upper power  bound of each generator
    P_1_text = P_i_min_all(1:D-1);
    P_2_text = P_i_max_all(1:D-1);
    for i=1:maxg             
        for j=1:sizepop      
            w= 0.4;                           %Inertia weight
            if j==1 
                if fitnesspbest(end)<fitnesspbest(2), lbest=pBest(end,:);
                else, lbest=pBest(2,:);
                end
            elseif j==sizepop
                 if fitnesspbest(sizepop-1)<fitnesspbest(1),lbest=pBest(sizepop-1);
                 else, lbest=pBest(1,:);
                 end
            else
                 if fitnesspbest(j-1)<fitnesspbest(j+1), lbest=pBest(j-1,:);
                 else, lbest=pBest(j+1,:);
                 end
            end       
            V(j,:)=w*V(j,:)+c1(i)*rand(1,D-1).*(pBest(j,:)-pop2(j,:))+c2(i)*rand(1,D-1).*(lbest-pop2(j,:));
%             V(j,:)=w*V(j,:)+c1(i)*rand*(pBest(j,:)-pop2(j,:))+c2(i)*rand*(gBest-pop2(j,:));
            V(j,find(V(j,:)<Vmin))=Vmin;
            V(j,find(V(j,:)>Vmax))=Vmax;
            pop2(j,:)=pop2(j,:)+V(j,:);
            ind_1 = pop2(j,:) < P_1_text;
            ind_peak = find(ind_1 == 1);
            pop2(j,ind_peak) = P_1_text(ind_peak);
            ind_2 = pop2(j,:) > P_2_text;
            ind_peak2 = find(ind_2 == 1);
            pop2(j,ind_peak2) = P_2_text(ind_peak2);
            pop2(find(pop2 < 0)) = 1;
            pop_exi = Total - sum(pop2(j,:));
            pop_sum = [pop2(j,:), pop_exi];
            if (pop_exi >  P_i_max_all(end)) || (pop_exi < P_i_min_all(end))
                fitness(j)= 1e10;
            else
                fitness(j)= power_allocation_6D(pop_sum);
            end
            if fitness(j) < fitnesspbest(j)
                pBest(j,:)= pop_sum(1:D-1);
                fitnesspbest(j)=fitness(j);
            end
            if fitness(j)<fitnessgbest
                gBest= pop_sum(1:D-1);
                fitnessgbest=fitness(j);
            end
        end
        result(i)=fitnessgbest;    
    end
    Y(iii,N_add+1:end) = result;
end
save('LPSO_EVOLER_6D','Y')
