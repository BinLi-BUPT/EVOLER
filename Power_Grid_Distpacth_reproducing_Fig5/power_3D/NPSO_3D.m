% clc;
clear all;
% close all;
c1 = 2;
c2 = 2;  
maxg = 1000;        
sizepop = 20;      
D= 3;               

Vmax = 20;
Vmin = -20;    
Total = 850;

P_i_min_all = [100,100,50];
P_i_max_all = [600,400,200];
popmax=P_i_max_all(1:D-1);
popmin=P_i_min_all(1:D-1);
DMin=P_i_min_all(end);
DMax=P_i_max_all(end);
P_1_text = P_i_min_all(1:D-1);
P_2_text = P_i_max_all(1:D-1);
pop = [];
pop2 = [];

for iii = 1 : 50
iii
for i=1:sizepop 
    while(1)
        pop2(i,:)= popmin(1:D-1)+(popmax(1:D-1)-popmin(1:D-1)).*rand(1,D-1);

        ind_1 = pop2(i,:) < P_1_text;
        ind_peak = find(ind_1 == 1);
        pop2(i,ind_peak) = P_1_text(ind_peak);

        ind_2 = pop2(i,:) > P_2_text;
        ind_peak2 = find(ind_2 == 1);
        pop2(i,ind_peak2) = P_2_text(ind_peak2);
        pop2(find(pop2 < 0)) = 1;
        if (Total - sum(pop2(i,:)))<=DMax && (Total - sum(pop2(i,:)))>=DMin, break;end
    end
    pop(i,:)= ([pop2(i,:), Total - sum(pop2(i,:))]);
    V(i,:)=Vmax*rands(1,D-1);            
    fitness(i)= quad_cf(pop(i,:));       
end
pBest= pop2;                            
fitnesspbest=fitness;                  
[fitnessgbest bestindex]=min(fitness);
gBest=(pop2(bestindex,:));


[fitnessgbest bestindex]=min(fitness);
[~, bestindex2]=max(fitness);
g_id = pop(bestindex2,:);
p_id = pop;
[fitnessgbest bestindex] = min(fitness);
gbest=pop(bestindex,:);
pbest=pop;
fitnesspbest = fitness;
fitnessgworst = g_id;
fitnesspworst = fitness;

result = zeros(1,maxg);

for i=1:maxg           
    for j=1:sizepop     
        w= 0.3;           
        V(j,:)=V(j,:)+c1*rand(1,D-1).*(pBest(j,:)-pop2(j,:))+c2*rand(1,D-1).*(gBest-pop2(j,:));
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
            fitness(j)= quad_cf(pop_sum); 
        end

%         pop_exi(find(pop_exi > P_i_max_all(end))) = P_i_max_all(end);
%         pop_exi(find(pop_exi < P_i_min_all(end))) = P_i_min_all(end);

        if fitness(j) < fitnesspbest(j)
            pBest(j,:)= pop_sum(1:D-1);
            fitnesspbest(j)=fitness(j);
        end

        if fitness(j)<fitnessgbest
            gBest= pop_sum(1:D-1);
            fitnessgbest=fitness(j);
        end
        
        if fitness(j) > fitnesspworst(j)
            p_id(j,:) = pop(j,:);
            fitnesspworst(j)=fitness(j);
        end
        if fitness(j) > fitnessgworst
            g_id =  pop(j,:);
            fitnessgworst = fitness(j);
        end
        
    end
    result(i)=fitnessgbest;   
end
g_Best_end = [gBest, (Total - sum(gBest))];
Y(iii,:) = result;
% hold on
% plot(result,'r')
% xlabel('Genration');
% ylabel('Fitness')
end
% Y_2 = mean(Y,1);
% figure
% plot(Y_2,'r')
% xlabel('Genration');
% ylabel('Fitness')

save('NPSO_3D','Y')

