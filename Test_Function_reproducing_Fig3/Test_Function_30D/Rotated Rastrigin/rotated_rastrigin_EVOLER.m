% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
MaxDT = 500;                                   % Maximum number of generations
c1 = linspace(2,1.8,MaxDT);                    % Cognition learning factor c1
c2 = linspace(1.8,1.0,MaxDT);                  % Social learning factor c2
w= 0.5;                                        % Inertia weight
D= 30;                                         % Dimension

sizepop = 100;                                 % Population size
popmax = 5.12;                                 % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;

s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
num_all = 500;                                 % The number of simulation trials
Fitness_all = [];
for iiii = 1 : num_all
    iiii
    M = rot_matrix(D);                         % Rot_matrix:randomly generate one rotated D-dimension matrix: Input: D--dimension; Output: randomly rotated matrix
    %% Step 1: Low-rank Representation Learning based on the hierarchy reconstruction strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    Discret_leng = 3;                          % Discrete size of the original problem space;
    Trail_times = 10;                          % The times of hierarchy reconstruction strategy, to reduce the problem space search range
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);
   
    % step (i): reduce the problem space search range by exploiting the hierarchy reconstruction strategy
    for ii  = 1 : Trail_times
        [~,pop_div1,pop_div2] = Search_30D_Initial(max_ite,min_ite,s-1,Discret_leng,M);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    
    % step (ii):identification of the global optimum by exploiting the hierarchy reconstruction strategy
    %------ Recon_30D_Initial: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    Trail_times2 = 8;                          % The times of hierarchy reconstruction strategy based on low-rank sampling
    for ii  = 1 : Trail_times2
        [~,pop_div1,pop_div2] = Recon_30D_Initial(max_ite,min_ite,s,Discret_leng,M);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end 
    [global_optimum] =  Recon_30D_Initial(max_ite,min_ite,s,Discret_leng+1,M);
    N_add = round(Trail_times*(5*Discret_leng^6*(s-1)^4-5*(s-1)^6)/sizepop)+...
        round(Trail_times2*((4*(Discret_leng)*s^4-3*s^5 + (Discret_leng)*s^5)*2+...
       (4*(Discret_leng)*s^5-3*s^5 + (Discret_leng)^2*s^4)*2+ (3*(Discret_leng)*s^5-2*s^5)*2)/sizepop)+...
       round(((4*(Discret_leng+1)*s^4-3*s^5 + (Discret_leng+1)*s^5)*2+...
       (4*(Discret_leng+1)*s^5-3*s^5 + (Discret_leng+1)^2*s^4)*2+ (3*(Discret_leng+1)*s^5-2*s^5)*2)/sizepop);    % The overall required samples to reconstruct an attention subspace, which is used in the final plot. 

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
        fitness(i)=rotated_rest(pop(i,:),M);
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
            fitness(j)=rotated_rest(pop(j,:),M);
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
save('Pro_data_rotrest','Fit_pro')


