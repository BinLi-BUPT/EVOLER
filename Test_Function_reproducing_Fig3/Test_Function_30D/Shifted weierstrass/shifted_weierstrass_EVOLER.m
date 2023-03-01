% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
MaxDT = 500;                                   % Maximum number of generations
c1 = linspace(2,0.5,MaxDT);                    % Cognition learning factor c1
c2 = linspace(1.5,2,MaxDT);                    % Social learning factor c2
w= 0.5;                                        % Inertia weight
D= 30;                                         % Dimension
sizepop = 100;                                 % Population size
popmax=0.5;                                    % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;


s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
num_all = 500;                                 % The number of simulation trials
Fitness_all = [];
for iiii = 1 : num_all
    iiii
    %% Step 1: Low-rank Representation Learning based on the hierarchy reconstruction strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    shift_op = (-0.1+0.2*rand(1,D));           % Shifted global optimum, each dimension range [-0.1,0.1]
    Trail_time = 28;                           % The times of hierarchy reconstruction
    Discret_leng =3;                           % Discrete size of the original problem space;
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);
    % identification of the global optimum by exploiting the hierarchy reconstruction strategy
    %------ Recon_30D_Initial_shift: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    for ii = 1 : Trail_time
        [global_optimum,pop_div1,pop_div2] = Recon_30D_Initial_shift(max_ite,min_ite,s,Discret_leng,D,shift_op);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    N_add = round(Trail_time*((4*Discret_leng*s^4-3*s^5 + Discret_leng^2*s^5)*2+...
        (4*Discret_leng*s^5-3*s^5 + Discret_leng^2*s^4)*2+ (3*Discret_leng*s^5-2*s^5) *2)/sizepop);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot
     
    %% Step 2: Evolutionary PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum.
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
        fitness(i)=shifted_weierstrass(pop(i,:),shift_op);
    end
    %----------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop; 
    fitnesspbest=fitness;
    %----------------------------------
    % Canonic PSO algorithm
    for i=1:MaxDT
        tic
        for j=1:sizepop
            V(j,:)=w*V(j,:)+c1(i)*rand*(pbest(j,:)-pop(j,:))+c2(i)*rand*(gbest-pop(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<-popmax))=-popmax;
            fitness(j)=shifted_weierstrass(pop(j,:),shift_op);
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
%     time_end(iiii) = time_1(iiii)+sum(time_2(1:900-N_add));
end

Fit_pro = Fitness_all';
Fit_pro = Fit_pro(:,1:900);
save('shifted_Pro_data_wei','Fit_pro')