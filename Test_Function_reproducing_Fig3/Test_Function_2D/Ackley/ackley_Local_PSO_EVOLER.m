% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
close all
clear all;
clc;
format long

%%
%-----------initialization parameter setting-------------------------
MaxDT=500;                                  % Maximum number of generations
c1 = linspace(2,0.5,MaxDT);                 % Cognition learning factor c1
c2 = linspace(1.5,2,MaxDT);                 % Social learning factor c2
w= 0.3;                                     % Inertia weight
D = 2;                                      % Dimension
sizepop=50;                                 % Population size

popmax=32;                                  % Variable range [popmin: popmax]
popmin=-32;                                 % Variable range [popmin: popmax]
Vmax=0.1*popmax;
Vmin=0.1*popmin;

Discret_leng = 80;                          % Discrete size of the original problem space; M=N=Discret_leng;
x_axis = linspace(popmin,popmax,Discret_leng);              
y_axis = linspace(popmin,popmax,Discret_leng);            

Discret_leng_now = 1:1:Discret_leng;
Downsample_factor = 4;
Discret_leng_int = Discret_leng_now(1:Downsample_factor:end); % Down-sampling with the factor 4

% The precomputed original problem space, for the subsequent sampling and reconstruction; 
% This is used to evaluate the reconstruction residual error of the learned low-rank representation. 
for iii=1:1:length(Discret_leng_int)
    for jjj=1:1:length(Discret_leng_int)
        Z(iii,jjj)=ackley([x_axis(Discret_leng_int(iii)),y_axis(Discret_leng_int(jjj))]);    % Ackley: a standard test function
    end
end

num_all = 500; % The number of simulation trials
for iiii = 1 : num_all  
    iiii
    %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    % step (i): Structured random sampling on the whole problem space,
    % i.e., sampling s rows and columns of the original discrete problem space;
    s = 16;                                            % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).  
    Index_c = randperm(length(Discret_leng_int),s);    % Random column index
    Index_r = randperm(length(Discret_leng_int),s);    % Random row index
    C1 = Z(:,Index_c);
    R1 = Z(Index_r,:);
    
    %----------------------------------
    % step (ii): Reconstruction of the approximate problem space, via the special form: \hat(Z) = C1 * U1 * R1;
    U0 = Z(Index_r,Index_c);
    [U_u,S_u,V_u]=svd(U0);
    Uu=U_u';
    s_0= s;
    U1 = V_u(:,1:s_0) * pinv(S_u(1:s_0,1:s_0)) * Uu(1:s_0,:); % Determining the central matrix U1.
    Z_est = C1*U1*R1;                                         % Reconstructing the problem space \hat(Z).
    
    %----------------------------------
    % step (iii): Identification of the global optimum in a representation space, and determination of the attention subspace     
    [X11,Y11] = meshgrid(Discret_leng_int);
    [X12,Y12] = meshgrid(Discret_leng_now);        
    Z_end = interp2(X11, Y11, Z_est, X12, Y12);               % Optional when the variable range is large, which helps to further reduce the number of samples.
    Z_end(Discret_leng_int,Discret_leng_int) = Z_est;         
    [m0,n0]=find(Z_end == min(min(Z_end)));                   % The center of the attention subspace, i.e., the global optimum in the representation space
    R = abs(x_axis(2) - x_axis(1));                           % The radius of the attention subspace
    
    N_add = round((s*length(Discret_leng_int)+s*length(Discret_leng_int)-s*s)/sizepop);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot. 
    
   %% Step 2: Evolutionary PSO method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
   %----------------------------------
   % Initialization of population, around the identified attention subspace
    x_po1 = x_axis(m0);
    y_po1 = y_axis(n0);
    see2(iiii,:) = [x_po1,y_po1];     
    theat_all = 2*pi*rand(sizepop,1);
    for i=1:sizepop
        if i == 1
            pop(i,:)= [x_axis(m0),y_axis(n0)];
        else
            pop(i,:)= [R*cos(theat_all(i))+x_axis(m0),R*sin(theat_all(i))+y_axis(n0)];
        end
        V(i,:)=Vmax*randn(1,D); %
        fitness(i)=ackley(pop(i,:));
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
            V(j,:)=w*V(j,:)+c1(i)*rand(1,D).*(pbest(j,:)-pop(j,:))+c2(i)*rand(1,D).*(lbest-pop(j,:));            V(j,find(V(j,:)>Vmax))=Vmax;
            V(j,find(V(j,:)<Vmin))=Vmin;
            pop(j,:)=pop(j,:)+V(j,:);
            pop(j,find(pop(j,:)>popmax))=popmax;
            pop(j,find(pop(j,:)<popmin))=popmin;
            fitness(j)=ackley(pop(j,:));
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
    Y = 100*ones(1,MaxDT);     % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy(1:end-N_add);
    
    Fitness_all(:,iiii)=Y;    
end

% figure
% semilogy(Fitness_all,'Color',[1 0.75 0.8])
% hold on
% plot(mean(Fitness_all'))
% xlabel('Genration');
% ylabel('Fitness')
des_val= 1e-15;  %If the output result is less than the threshold 'des_val', the global optimum is found  
Fit_pro = Fitness_all';
save('Pro_Local_data_ack','Fit_pro')

figure
semilogy(mean(Fitness_all'),'LineWidth',2);
xlabel('Genration');
ylabel('Fitness')
legend('Local PSO EVOLER');


da_1_Pro = load('Pro_Local_data');  % run the file 'rosen_EVOLER' to generate this data
ff_end15 = da_1_Pro.Fit_pro(1:num_all,:);
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
set(h15,'MarkerFaceColor',get(h15,'color'));
% legend('Local PSO EVOLER');
