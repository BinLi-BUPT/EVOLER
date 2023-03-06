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
w= 0.3;                                        % Inertia weight
D= 2;                                          % Dimension

sizepop = 50;                                  % Population size
popmax = 5;                                    % Variable Range [-popmax: popmax]
Vmax=0.1*popmax;
Vmin=-0.1*popmax;
Discret_leng = 40;                            % Discrete size of the original problem space;  
x_axis = linspace(-popmax,popmax,Discret_leng);  
y_axis = linspace(-popmax,popmax,Discret_leng);  

Discret_leng_now = 1:1:Discret_leng; 
Downsample_factor = 2;                         % Down-sampling with the factor 
Discret_leng_int = Discret_leng_now(1:Downsample_factor:end); 
for iii=1:1:length(Discret_leng_int)
    for jjj=1:1:length(Discret_leng_int)
        Z(iii,jjj)=hybrid_func2([x_axis(Discret_leng_int(iii)),y_axis(Discret_leng_int(jjj))]);  % hybrid_func2: a standard test function
    end
end
num_all = 500;                                 % The number of simulation trials
Fitness_all = [];
for iiii = 1 : num_all
    iiii    
    %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    % step (i): Structured random sampling on the whole problem space,
    % i.e., sampling s rows and columns of the original discrete problem space;
    s = 4;                                     % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m). 
    Index_c = randperm(length(Discret_leng_int),s);        % Randomly column index set
    Index_r = randperm(length(Discret_leng_int),s);        % Randomly row index set
    C1 = Z(:,Index_c);                           
    R1 = Z(Index_r,:);                             

    %----------------------------------
    % step (ii): Reconstruction of the approximate problem space, via the special form: \hat(Z) = C1 * U1 * R1;
    U0 = Z(Index_r,Index_c);
    [U_u,S_u,V_u]=svd(U0);
    Uu=U_u';
    s_0= s-1;
    U1 = V_u(:,1:s_0) * pinv(S_u(1:s_0,1:s_0)) * Uu(1:s_0,:); % Determining the central matrix U1.
    Z_est = C1*U1*R1;                                         % Reconstructing the problem space \hat(Z).
    
    %----------------------------------
    % step (iii): Identification of the global optimum in a representation space, and determination of the attention subspace     
    [X11,Y11] = meshgrid(Discret_leng_int);
    [X12,Y12] = meshgrid(Discret_leng_now);        
    Z_end = interp2(X11, Y11, Z_est, X12, Y12);               % Optional when the variable range is very large, which helps to further reduce the number of samples.
    Z_end(Discret_leng_int,Discret_leng_int) = Z_est;          
    [m0,n0]=find(Z_end == min(min(Z_end)));                  
    global_optimum = [x_axis(m0(1)),y_axis(n0(1))];           % The center of the attention subspace, i.e., the global optimum in the representation space
    N_add = round((2*s*length(Discret_leng_int)-s*s)/sizepop);% The overall required samples to reconstruct an attention subspace, which is used in the final plot.
    R = abs(x_axis(2) - x_axis(1));                           % The radius of the attention subspace
    theat_all = 2*pi*rand(sizepop,1);
    %% Step 2: Evolutionary PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
    %----------------------------------
    % Initialization of population, around the identified attention subspace
    for i=1:sizepop 
        if i == 1
            pop(i,:)= global_optimum;
        else
            pop(i,:)= [R*cos(theat_all(i))+global_optimum(1),R*sin(theat_all(i))+global_optimum(2)];
        end
        V(i,:) = rand(1,D); % 
        fitness(i)=hybrid_func2(pop(i,:));
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
            fitness(j)=hybrid_func2(pop(j,:));
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
save('Pro_data','Fit_pro')


