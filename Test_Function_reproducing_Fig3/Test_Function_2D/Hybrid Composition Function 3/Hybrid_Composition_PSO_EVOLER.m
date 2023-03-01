% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;
close all

%%
%-----------initialization parameter setting-------------------------
MaxDT = 200;                                   % Maximum number of generations
c1 = linspace(2,0.5,MaxDT);                    % Cognition learning factor c1
c2 = linspace(1.5,2,MaxDT);                    % Social learning factor c2
w= 0.3;                                        % Inertia weight
D= 2;                                          % Dimension
sizepop = 50;                                  % Population size
popmax = 5;                                    % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;
num_all =500;                                  % The number of simulation trials
Fitness_all = [];

for iiii = 1 : num_all
    iiii    
    %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    Discret_leng = 20;                         % Discrete size of the original problem space;
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);
    
    % step (i): reduce the problem space search range based on full sampling strtegy
    [~,pop_div1,pop_div2] = Search_2D_Initial(min_ite,max_ite,Discret_leng);
    max_ite = pop_div1;
    min_ite = pop_div2;

    % step (ii): identification of the global optimum based on low-rank sampling strtegy
    %------ Recon_30D_Initial: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    Trail_time = 3;                           % The repeating times of dichotomy strategy based on low-rank sampling
    Discret_leng1 = 3;
    s = 2;                                    % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
    for ii = 1 : Trail_time
        [global_optimum,pop_div1,pop_div2] = Recon_2D_Initial(min_ite,max_ite,s,Discret_leng1);   % global_optimum is the center of the attention subspace, i.e., the global optimum in the representation space
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    R = abs(pop_div2(1) - pop_div1(1));       % The radius of the attention subspace  
    N_add = round((Discret_leng*Discret_leng + Trail_time*(Discret_leng1*s*2-s^2))/sizepop);   % The overall required samples to reconstruct an attention subspace, which is used in the final plot.
    
    %% Step 2: Evolutionary PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
    %----------------------------------
    % Initialization of population, around the identified attention subspace
    theat_all = 2*pi*rand(sizepop,1);
    for i=1:sizepop 
        if i == 1
            pop(i,:)= global_optimum;
        else
            pop(i,:)= [R*cos(theat_all(i))+global_optimum(1),R*sin(theat_all(i))+global_optimum(2)];
        end
        V(i,:) = rand(1,D); % 
        fitness(i)=hybrid_func3(pop(i,:));
    end
    %----------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop;
    fitnesspbest=fitness;
    
    %----------------------------------
    % Canonic PSO algorithm
    for i=1:MaxDT
        for j=1:sizepop
            V(j,:)=w*V(j,:)+c1(i)*rand*(pbest(j,:)-pop(j,:))+c2(i)*rand*(gbest-pop(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<-popmax))=-popmax;
            fitness(j)=hybrid_func3(pop(j,:));
            if fitness(j)<fitnesspbest(j)
                pbest(j,:)=pop(j,:);
                fitnesspbest(j)=fitness(j);
            end
            if fitness(j)<fitnessgbest
                gbest=pop(j,:);
                fitnessgbest=fitness(j);
            end
        end
        yy(i)= fitnessgbest;   
    end
    
    Y = 100*ones(1,N_add+MaxDT);        % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy;
    Fitness_all(:,iiii)=Y;
end
figure
semilogy(Fitness_all,'Color',[1 0.75 0.8])
hold on
plot(mean(Fitness_all'))
xlabel('Genration');
ylabel('Fitness')
 
Fit_pro = Fitness_all';
save('Pro_data_hybrid_3','Fit_pro')


