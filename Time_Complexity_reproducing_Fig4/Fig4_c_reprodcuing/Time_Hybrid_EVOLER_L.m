% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
des_val= 8e-4;                                 % If the output result is less than the threshold 'des_val', the global optimum is found  
maxGen = 500;                                  % Maximum number of generations
c1 = linspace(2,0.5,maxGen);                   % Cognition learning factor c1
c2 = linspace(1.5,2,maxGen);                   % Social learning factor c2
w= 0.5;                                        % Inertia weight
D= 30;                                         % Dimension

sizepop = 100;                                 % Population size
popmax = 5;                                    % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;

s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m). 
num_all = 10;                                 % The number of simulation trials
Fitness_all = [];
for iiii = 1 : num_all
    iiii    
    %% Step 1: Low-rank Representation Learning based on the dichotomy strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    Discret_leng = 3;                          % Discrete size of the original problem space;                              
    Trail_times = 18;                          % The repeating times of dichotomy strategy based on full sampling
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);
    PDF_cal = 0;
    tic
    % step (i): reduce the problem space search range based on full sampling strtegy
    for ii  = 1 : Trail_times
        [~,pop_div1,pop_div2] = Search_30D_Initial(max_ite,min_ite,s-1,Discret_leng);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    
    % step (ii): identification of the global optimum based on low-rank sampling strtegy
    Discret_leng1 = Discret_leng+1;            % Discrete size of the original problem space for low-rank sampling       
    Trail_times2 = 3;                          % The repeating times of dichotomy strategy based on low-rank sampling
    %------ Recon_30D_Initial: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    for ii = 1 : Trail_times2
        [global_optimum,pop_div1,pop_div2] = Recon_30D_Initial(max_ite,min_ite,s,Discret_leng1);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end    
    
    N_add = round(Trail_times*(5*Discret_leng^6*(s-1)^4-5*(s-1)^6)/sizepop)+...
       round(Trail_times2*((4*Discret_leng1*s^5-3*s^5 + Discret_leng1^2*s^4)*2+ (3*Discret_leng1*s^5-2*s^5) *2)/sizepop);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot

   %% Step 2: Evolutionary Local PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
    %----------------------------------
    % Initialization of population, around the identified attention subspace
    R=(pop_div1-pop_div2);
    for i=1:sizepop 
        if i == 1
            pop(i,:)= global_optimum;
        else
            pop(i,:)= global_optimum +R.*rand(1, D);
        end
        V(i,:) = rand(1,D); % 
        fitness(i)=hybrid_func1(pop(i,:));
    end
    %----------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop;
    fitnesspbest=fitness;
    time_1(iiii) = toc;
    %----------------------------------
    % Canonic Local PSO algorithm
    for i=1:maxGen
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
            fitness(j)=hybrid_func1(pop(j,:));
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
        if fitnessgbest > des_val
            PDF_cal = PDF_cal + 1;
        end
        time_2(i) = toc;
    end
    Y = 100*ones(1,N_add+maxGen);        % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy;
    Fitness_all(:,iiii)=Y;
    time13(iiii) = time_1(iiii)+sum(time_2(1:900-N_add));
end
save time_evoler_hybrid_L time13