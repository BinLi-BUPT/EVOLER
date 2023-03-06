% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
des_val= 1e-15;                                % If the output result is less than the threshold 'des_val', the global optimum is found  
MaxDT = 600;                                   % Maximum number of generations
c1 = linspace(2,0.5,MaxDT);                    % Cognition learning factor c1
c2 = linspace(1.5,2,MaxDT);                    % Social learning factor c2
w= 0.5;                                        % Inertia weight
D= 30;                                         % Dimension
sizepop = 100;                                 % Population size
popmax= 20;                                    % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;

s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
num_all = 101;                                 % The number of simulation trials
Fitness_all = [];
for iiii = 1 : num_all
    iiii
    %% Step 1: Low-rank Representation Learning based on the dichotomy strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    shift_op = (-2+4*rand(1,D));               % Shifted global optimum, range [-2,2]
    Trail_time = 19;                           % The repeating times of dichotomy strategy based on low-rank sampling
    Discret_leng = 4;                          % Discrete size of the original problem space;      
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);
    
    % identification of the global optimum based on low-rank sampling strtegy
    %------ Recon_30D_Initial_shift: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    tic
    for ii = 1 : Trail_time
        [global_optimum,pop_div1,pop_div2] = Recon_30D_Initial_shift(max_ite,min_ite,s,Discret_leng,D,shift_op);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    N_add = round(Trail_time*((4*Discret_leng*s^4-3*s^5 + Discret_leng^2*s^5)*2+...
        (4*Discret_leng*s^5-3*s^5 + Discret_leng^2*s^4)*2+ (3*Discret_leng*s^5-2*s^5) *2)/sizepop);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot
     
    %% Step 2: Evolutionary Local PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum.
    %----------------------------------
    % Initialization of population, around the identified attention subspace
    R=(pop_div1-pop_div2);                     % The radius of the attention subspace  
    for i=1:sizepop 
        if i == 1
            pop(i,:)= global_optimum;
        else
            pop(i,:)= global_optimum +R.*rand(1, D);
        end
        V(i,:) = rand(1,D); %
        fitness(i)=shifted_levy(pop(i,:),shift_op);
    end
    %----------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop; 
    fitnesspbest=fitness;
    time_1(iiii) = toc;
    %----------------------------------
    % Canonic Local PSO algorithm
    for i=1:MaxDT
        tic
        for j=1:sizepop
            if j==1 
                if fitnesspbest(end)<fitnesspbest(2), lbest=pbest(end,:);
                else, lbest=pbest(2,:);
                end
            elseif j==sizepop
                 if fitnesspbest(sizepop-1)<fitnesspbest(1),lbest=pbest(sizepop-1);
                 else, lbest=pbest(1,:);
                 end
            else
                 if fitnesspbest(j-1)<fitnesspbest(j+1), lbest=pbest(j-1,:);
                 else, lbest=pbest(j+1,:);
                 end
            end
            V(j,:)=w*V(j,:)+c1(i)*rand(1,D).*(pbest(j,:)-pop(j,:))+c2(i)*rand(1,D).*(lbest-pop(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<-popmax))=-popmax;
            fitness(j)=shifted_levy(pop(j,:),shift_op);
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
        time_2(i) = toc;
    end
    Y = 100*ones(1,MaxDT+N_add);        % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy;
    Fitness_all(iiii,:)=Y;
    time_evoler_l(iiii) = time_1(iiii)+sum(time_2(1:1200-N_add));
end
% time_evoler_l = time_evoler_l(2:end);
% Fit_pro = Fitness_all';
Fit_pro = Fitness_all(2:end,1:900);
save('Local_data_levy','Fit_pro')
% save time_evoler_l time_2

