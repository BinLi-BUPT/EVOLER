%-------------Traditional PSO algorithm----------%
close all
clear all;
clc;
format long;

c1=1.8;               % Cognition learning factor c1   
c2=1.8;               % Cognition learning factor c1          
w= 0.3;               % Inertia weight
D= 2;                 % Dimension       
MaxDT=500;            % population size
N=50;                  
Vmax=2;
Vmin=-2;
popmax=500;
popmin=-500;
Num_exper=500;
Fitness_all = zeros(MaxDT,Num_exper);

for kkkkk=1:Num_exper
        
    for i=1:N
        pop(i,:)=popmin+(popmax-popmin)*rand(1,2); 
        V(i,:)=rand(1,2);
        fitness(i)=schewfel(pop(i,:));
    end    
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop;
    fitnesspbest=fitness;
    
    for i=1:MaxDT
        for j=1:N
            V(j,:)=w*V(j,:)+c1*rand*(pbest(j,:)-pop(j,:))+c2*rand*(gbest-pop(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<popmin))=popmin;
            fitness(j)=schewfel(pop(j,:));
            if fitness(j)<fitnesspbest(j)
                pbest(j,:)=pop(j,:);
                fitnesspbest(j)=fitness(j);
            end
            
            if fitness(j)<fitnessgbest
                gbest=pop(j,:);
                fitnessgbest=fitness(j);
            end
            
        end
        yy(i)=fitnessgbest;
        
    end
    Fitness_all(:,kkkkk)=yy;
    kkkkk
end

%---Fig.2-d: the convergence of classical evolutionary computing PSO with 50 particles. 
figure(111)
semilogy(Fitness_all,'Color',[0.7,0.7,1.00])
hold on
semilogy(mean(Fitness_all'))
xlabel('Genration');
ylabel('Fitness')

%---Fig. 2-(f): the convergence Generations of traditional PSO method
Convergence_point = [];
for kkkkk=1:Num_exper
    fitness_current = Fitness_all(:,kkkkk);
    x_convergent = find(fitness_current<5e-5);
    if (~isempty(x_convergent))
        Convergence_point = [Convergence_point; min(x_convergent)];
    end
end
[Yy1,Xx1] = hist(Convergence_point,20);
figure
hold on
bar(Xx1,Yy1)
xlabel('Convergence Generations');
ylabel('Histogram')
% figure
% subplot(2,1,1)
% xxx=Fitness_all(100,:);
% [xx1,yy1]=hist(xxx,200);
% plot(yy1,xx1)
% xlabel('Fitness');
% ylabel('Histogram');
% 
% subplot(2,1,2)
% xxx=Fitness_all(500,:);
% [xx1,yy1]=hist(xxx,200);
% plot(yy1,xx1)
% xlabel('Fitness');
% ylabel('Histogram');


 