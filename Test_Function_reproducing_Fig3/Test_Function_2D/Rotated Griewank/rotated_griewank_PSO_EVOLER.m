% Copyright 2022, All Rights Reserved
% Code by Bin Li, Ziping Wei
clear all;
clc;
format long;
close all

%%
%-----------initialization parameter setting-------------------------
MaxDT=500;                                      % Maximum number of generations
c1 = linspace(2,0.5,MaxDT);                     % Cognition learning factor c1
c2 = linspace(1.5,2,MaxDT);                     % Social learning factor c2
w= 0.3;                                         % Inertia weight
D=2;                                            % Dimension
N=50;                                           % Population size

popmax= 600;                                    % Variable Range [popmin: popmax]
popmin=-popmax;                                 % Variable Range [popmin: popmax]
Vmax=0.1*popmax;
Vmin=0.1*popmin;

Discret_leng = 4000;                            % Discrete size of the original problem space;
x_axis = linspace(popmin,popmax,Discret_leng);  
y_axis = linspace(popmin,popmax,Discret_leng);  

Discret_leng_now = 1:1:Discret_leng;
Downsample_factor = 20;
Discret_leng_int = Discret_leng_now(1:Downsample_factor:end);  % Down-sampling with the factor 80
Discret_leng_int_num = length(Discret_leng_int);

num_all = 500; % The number of simulation trials
for iiii = 1 : num_all
    iiii
    % The precomputed original problem space, for the subsequent sampling and reconstruction;
    % This is used to evaluate the reconstruction residual error of the learned low-rank representation.
    M = rot_matrix(D); % rot_matrix:randomly generate one rotated D-dimension matrix: Input: D--dimension; Output: randomly rotated matrix
    for iii=1:1:length(Discret_leng_int)
        for jjj=1:1:length(Discret_leng_int)
            Z(iii,jjj)=rotated_griewank([x_axis(Discret_leng_int(iii)),y_axis(Discret_leng_int(jjj))],M); % rotated_griewank: a standard test function
        end   
    end
    
   %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    % step (i): Structured random sampling on the whole problem space,
    % i.e., sampling s rows and columns of the original discrete problem space;
    s = 7;                                          % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m). 
    Index_c = randperm(length(Discret_leng_int),s); % Randomly column index set
    Index_r = randperm(length(Discret_leng_int),s); % Randomly row index set
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
    Z_end = interp2(X11, Y11, Z_est, X12, Y12);        % Optional when the variable range is very large, which helps to further reduce the number of samples.
    Z_end(Discret_leng_int,Discret_leng_int) = Z_est;  
    [m0,n0]=find(Z_end == min(min(Z_end)));            % The center of the attention subspace, i.e., the global optimum in the representation space
    R = abs(x_axis(2) - x_axis(1));                    % The radius of the attention subspace
    
    N_add = round((s*(Discret_leng_int_num)+s*(Discret_leng_int_num)-s*s)/N);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot.  
    
   %% Step 2: Evolutionary PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
   %----------------------------------
   % Initialization of population, around the identified attention subspace
    m=m0(1);
    n=n0(1);
    see(iiii,:) = [m,n];
    theat_all = 2*pi*rand(N,1);
    for i=1:N
        if i == 1
            pop(i,:)= [x_axis(m), y_axis(n)];
        else
            pop(i,:)= [R*cos(theat_all(i))+x_axis(m0),R*sin(theat_all(i))+y_axis(n0)];
        end
        V(i,:)= rand(1,D); %
        fitness(i)=rotated_griewank(pop(i,:),M);
    end
    %-------------------------------------
    [fitnessgbest bestindex]=min(fitness);
    gbest=pop(bestindex,:);
    pbest=pop;
    fitnesspbest=fitness;
    
    %----------------------------------
    % Canonic PSO algorithm
    for i=1:MaxDT
        for j=1:N
            V(j,:)=w*V(j,:)+c1(i)*rand*(pbest(j,:)-pop(j,:))+c2(i)*rand*(gbest-pop(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<popmin))=popmin;
            fitness(j)=rotated_griewank(pop(j,:),M);
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
    Y = 100*ones(1,MaxDT);    % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy(1:end-N_add);
    Fitness_all(:,iiii)=Y;
    
end

% figure
% semilogy(Fitness_all,'Color',[1 0.75 0.8])
% hold on
% plot(mean(Fitness_all'))
% xlabel('Genration');
% ylabel('Fitness')
Fit_pro = Fitness_all';
save('Pro_data_rotgriewank','Fit_pro')
