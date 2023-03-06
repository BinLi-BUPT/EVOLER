% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
close all;

%%
%-----------initialization parameter setting-------------------------
c1 = 2;                                  % Cognition learning factor c1
maxg = 400;                              % maximum number of generations
c2 = linspace(1.8,1.9,maxg);             % Social learning factor c2
sizepop = 20;                            % population size
D= 2;                

Vmax = 20;
Vmin = -20;
P_i_min = [100,100,50];                  % lower power limitation of each generator
P_i_max = [600,400,200];                 % upper power limitation of each generator
popmax = 500;                            % variable Range [popmin: popmax]
popmin = 50;                             % variable Range [popmin: popmax]

Discret_leng = 10;                       % discrete size of the original problem space;
x_axis = linspace(popmin,popmax,Discret_leng);
x_axis = x_axis(2:end);
y_axis = linspace(popmin,popmax,Discret_leng);
y_axis = y_axis(2:end);
z_axis = linspace(popmin,popmax,Discret_leng);
z_axis = z_axis(2:end);
Discret_leng2 = length(z_axis);
Total = 850;                             % total power limitation

% The precomputed original problem space, for the subsequent sampling and reconstruction; 
% This is used to evaluate the reconstruction residual error of the learned low-rank representation. 
for iii=1:1:length(x_axis)
    for jjj=1:1:length(y_axis)
        for kkk=1:1:length(z_axis)
            x1 = x_axis(iii);
            x2 = y_axis(jjj);
            x3 = z_axis(kkk);
            % function quad_cf: to calculate the cost of D=3 generators;
            % Input: the power of each generator; Output: cost;
            % see https://alroomi.org/economic-dispatch.
            Z(iii,jjj,kkk)=quad_cf([x1,x2,x3]);
        end
    end
end

% Singular values of the problem space (D = 3): Reproducing Figure 4-b
s11 = tenmat(Z,1);
[~,sig_1,~] = svd(double(s11),'econ');
figure
semilogy(diag(sig_1))
s11 = tenmat(Z,2);
[~,sig_1,~] = svd(double(s11),'econ');
hold on
semilogy(diag(sig_1))
s11 = tenmat(Z,3);
[~,sig_1,~] = svd(double(s11),'econ');
hold on
semilogy(diag(sig_1))

Z_ten = tensor(Z);
A_all= [];
N_all = 50;                  % The number of simulation trials
c = 2;                       % the sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).  
N_add = round((length(x_axis)*(c)*3 - 2*(c)*(c)*(c))/sizepop);% the overall required samples to reconstruct an attention subspace, which is used in the final plot.
Y = 8.5e3*ones(N_all,N_add+maxg);
for sss = 1 : N_all
    
   %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    % step (i): Structured random sampling on the whole problem space;
    A = Z;
    n_all = size(A);
    n1 = n_all(1);
    n2 = n_all(2);
    n3 = n_all(3);
    c_ll = [randperm(n1,c)];      % sampling index set of 1 dimension
    c_l2 = [randperm(n2,c)];      % sampling index set of 2 dimension
    c_l3 = [randperm(n3,c)];      % sampling index set of 3 dimension
    R =tensor(A(c_ll,c_l2,c_l3)); % sampling core tensor 
    C1_pr = A(:,c_l2,c_l3);       % sampling first small tensor 
    C2_pr = A(c_ll,:,c_l3);       % sampling second small tensor 
    C3_pr = A(c_ll,c_l2,:);       % sampling third small tensor 
    
    % step (ii): Reconstruction of the approximate high-dimensional problem space
    %--------------Tensor CUR recovery algorithm----------------------%
    % Input: C1_pr:sampling first small tensor
    % C2_pr:sampling second small tensor 
    % C3_pr:sampling third small tensor 
    % R: core tensor
    % c_ll ;  sampling index set of first dimension
    % c_l2  ; sampling index set of second dimension
    % c_l3 ;  sampling index set of third dimension
    % Output: A_result: the recovered representation space
    A_result = Tensor3_CUR(C1_pr,C2_pr,C3_pr,R,c_ll,c_l2,c_l3); 
    
    % step (iii): Identification of the global optimum in the representation space, and determination of the attention subspace 
    A_result2 = double(A_result);
    for iii=1:1:Discret_leng2
        for jjj=1:1:Discret_leng2
            for kkk=1:1:Discret_leng2
                x1 = x_axis(iii);
                x2 = y_axis(jjj);
                x3 = Total - x1 - x2;
                x3_o = z_axis(kkk);
                if x3_o == x3
                else
                    A_result2(iii,jjj,kkk) = 1e5;
                end
            end
        end
    end
    [m_all22] = find(tensor(A_result2) == collapse(tensor(A_result2),[1,2,3],@min));  % The center of the attention subspace, i.e., the global optimum in the representation space
    x_min2 = x_axis(m_all22(1));
    y_min2 = y_axis(m_all22(2));
    z_min2 = Total - x_min2 - y_min2;
    A_all(sss,:) = [x_min2,y_min2,z_min2]; 
    
    
   %% Step 2: Evolutionary PSO method, which is used to exploit the identified attention subspace and finally gain the global optimum. 
   % Initialization  of population, around the identified attention subspace
   pop = [];
    pop2 = [];
    for i=1:sizepop 
        if i == 1
            pop(i,:)= ([x_min2, y_min2, z_min2]);
            pop2(i,:)= ([x_min2, y_min2]);
        else
            x_pr = (rand*2)+x_min2;
            y_pr = (rand*2)+y_min2;
            pop(i,:)= ([x_pr, y_pr, Total - x_pr - y_pr]); 
            pop2(i,:) = ([x_pr, y_pr]);
        end
        V(i,:)=Vmax*rands(1,D);            
        fitness(i)= quad_cf([pop(i,:)]);      
    end
    pBest= pop2;                             
    fitnesspbest=fitness;                  
    [fitnessgbest bestindex]=min(fitness);
    gBest=(pop2(bestindex,:));
    result = zeros(1,maxg);
    for i=1:maxg             
        for j=1:sizepop      
            w= 0.3;           
            V(j,:)=w*V(j,:)+c1*rand*(pBest(j,:)-pop2(j,:))+c2(i)*rand*(gBest-pop2(j,:));
            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop2(j,:)=pop2(j,:)+V(j,:);
            % Limit the searching range 
            pop2(find(pop2(j,1) < P_i_min(1)),1) = P_i_min(1);
            pop2(find(pop2(j,2) < P_i_min(2)),2) = P_i_min(2);
            pop2(find(pop2(j,1) > P_i_max(1)),1) = P_i_max(1);
            pop2(find(pop2(j,2) > P_i_max(2)),2) = P_i_max(2);
            pop2(find(pop2 < 0)) = 1;
            
            pop_exi = Total - pop2(j,1) - pop2(j,2);
            pop_exi(find(pop_exi > P_i_max(3))) = P_i_max(3);
            pop_exi(find(pop_exi < P_i_min(3))) = P_i_min(3);
            pop_sum = [pop2(j,:), pop_exi];
            % calculating the fitness
            fitness(j)= quad_cf(pop_sum);
            % Update the optimal location of individual one
            if fitness(j) < fitnesspbest(j)
                pBest(j,:)= pop_sum(1:2);
                fitnesspbest(j)=fitness(j);
            end
            % Update the optimal location of group
            if fitness(j)<fitnessgbest
                gBest= pop_sum(1:2);
                fitnessgbest=fitness(j);
            end
        end
        result(i)=fitnessgbest; 
    end
    Y(sss,N_add+1:end) = result;
end
% figure
% plot([1:size(Y,2)],Y')
% xlabel('Genration')
% ylabel('Fitness')
% Y22 = mean(Y,1);
% hold on
% plot([1:size(Y,2)],Y22,'r')
save('EVOLER_3D','Y')