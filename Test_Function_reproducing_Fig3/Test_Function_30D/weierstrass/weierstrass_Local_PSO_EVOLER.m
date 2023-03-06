% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
des_val= 2e-14;                                % If the output result is less than the threshold 'des_val', the global optimum is found  
MaxDT = 200;                                   % Maximum number of generations
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
    Discret_leng1 = 3;                         % Discrete size of the original problem space;
    Trail_time1 = 49;                          % The times of hierarchy reconstruction
    max_ite = repmat(popmax,1,D);
    min_ite = repmat(-popmax,1,D);  
    
    % identification of the global optimum by exploiting the hierarchy reconstruction strategy
    %------ Recon_30D_Initial: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    for ii = 1 : Trail_time1
        [global_optimum,pop_div1,pop_div2] = Recon_30D_Initial(max_ite,min_ite,s,Discret_leng1);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    N_add =  round(Trail_time1*((4*Discret_leng1*s^4-3*s^5 + Discret_leng1^2*s^5)*2+...
        (4*Discret_leng1*s^5-3*s^5 + Discret_leng1^2*s^4)*2+ (3*Discret_leng1*s^5-2*s^5) *2)/sizepop); % The overall required samples to reconstruct an attention subspace, which is used in the final plot
     
    %% Step 2: Evolutionary Local PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum.
    %----------------------------------
    % Initialization of population, around the identified attention subspace
    R=(pop_div1-pop_div2);
    for i=1:sizepop
        if i == 1
            pop(i,:)= global_optimum;
        else
            pop(i,:)= global_optimum + R.*rand(1, D);
        end
        V(i,:) = rand(1,D); %
        fitness(i)=weierstrass(pop(i,:));
    end
    %----------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop; 
    fitnesspbest=fitness;
    
    %----------------------------------
    % Canonic Local PSO algorithm
    for i=1:MaxDT
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
            fitness(j)=weierstrass(pop(j,:));
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

Fit_pro = Fitness_all';
figure
semilogy(mean(Fitness_all'),'LineWidth',2);
xlabel('Genration');
ylabel('Fitness')
legend('Local PSO EVOLER');
ff_end15 = Fit_pro(1:num_all,:);
pro_15 = mean(ff_end15);
hold on
PSO_15 = [];
for ii = 1 :size(ff_end15,1)
    PDF_cal15 = 0;
    dp_now = ff_end15(ii,:);
    for jj = 1 : length(dp_now)
        if dp_now(jj) > des_val
            PDF_cal15 = PDF_cal15 + 1;
        end
    end
    PDF_num_LPSO(ii) = PDF_cal15;
    if ff_end15(ii,end-3) < des_val
        PSO_15(ii) = 1;
    end
end
X1 = PDF_num_LPSO;
[~,id_15] = find(PSO_15 == 1);
PDF_num_LPSO_15 = PDF_num_LPSO(id_15);
[y_pso,x_pso]=hist(PDF_num_LPSO_15,20);
if isempty(PDF_num_LPSO_15)
    PDF_num_LPSO_15 = 500;
    p15 = 0;
else
    p15 = length(PDF_num_LPSO_15)/num_all;
end
mean_various=mean(PDF_num_LPSO_15);
Probablity=p15;
figure
h15 = plot(mean_various,Probablity,'x','MarkerSize',12,'lineWidth',2,'color',[0.15,0.15,0.15]);
legend('Local PSO EVOLER');
