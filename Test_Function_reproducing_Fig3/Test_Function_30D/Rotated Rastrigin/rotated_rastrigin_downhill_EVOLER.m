% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
des_val= 1e-5;                                 % If the output result is less than the threshold 'des_val', the global optimum is found  
maxGen = 500;                                  % Maximum number of generations
D= 30;                                         % Dimension
sizepop = 100;                                 % Population size
popMax = 5.12;                                 % Variable Range [-popmax: popmax]
popMin=-popMax;

s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
num_all = 500;                                 % The number of simulation trials
Fitness_all = [];
proposed_downhill_16=zeros(1,num_all);
for iiii = 1 : num_all
    iiii
    M = rot_matrix(D);                         % Rot_matrix:randomly generate one rotated D-dimension matrix: Input: D--dimension; Output: randomly rotated matrix
    %% Step 1: Low-rank Representation Learning based on the hierarchy reconstruction strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    Discret_leng = 3;                          % Discrete size of the original problem space;
    Trail_times = 10;                          % The times of hierarchy reconstruction strategy, to reduce the problem space search range
    max_ite = repmat(popMax,1,D);
    min_ite = repmat(-popMax,1,D);
   
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

    %% Step 2: Downhill method, which is used to exploit the identified Attention Subspace and finally gain the global optimum.
    %----------------------------------
    % parameter setting
    sigma=0.5; % shrink coefficient
    alpha=1; % reflection coefficient
    beta=2; % expansion coefficient
    gamma=0.5;% contraction coefficient
    xMax=popMax;
    xMin=popMin;
    R=pop_div1 - pop_div2;
    theat_all = 3*pi*rand(D+1,1);
    % Build the initial simplex 
    for i=1:D+1
        if i==1
            S(:,i)=global_optimum;
        else
            S(:,i)=S(:,i-1);
            S(i-1,i)=S(i-1,i-1)+R(i-1);
        end
        S(find(S(:,i)<xMin),i)=xMin;
        S(find(S(:,i)>xMax),i)=xMax;
        fit_S(i)=rotated_rest(S(:,i)',M);
    end
    % Sort in descending order of fitness
    [fit_S,ind]=sort(fit_S,'descend');
    S_tmp=S;
    S=S_tmp(:,ind);
    PDF_cal16=N_add;
    yy16(1)=min(fit_S);
   for i=1:maxGen
        Vh=S(:,1);
        fit_Vh=fit_S(1);
        Vl=S(:,end);
        fit_Vl=fit_S(end);   
        Vx=mean(S(:,2:end),2);% centroid
        Vr=Vx+alpha.*(Vx-Vh);% reflection
        for d=1:D
            Vr(d)=min(xMax,Vr(d));
            Vr(d)=max(xMin,Vr(d));
        end
        fit_Vr=rotated_rest(Vr',M);%***********************
        if fit_Vr<=fit_Vl % Expansion:he reflection perform berter than the best position. 
            Ve=Vx+beta.*(Vr-Vx);
            for d=1:D
                Ve(d)=min(xMax,Ve(d));
                Ve(d)=max(xMin,Ve(d));
            end
            fit_Ve=rotated_rest(Ve',M);%******************
            if fit_Ve<fit_Vl
                S(:,1)=[];      S=[S,Ve];
                fit_S(1)=[];    fit_S=[fit_S,fit_Ve];
            else
                S(:,1)=[];      S=[S,Vr];
                fit_S(1)=[];    fit_S=[fit_S,fit_Vr];
            end
        elseif fit_Vr<fit_S(2) % The reflection performs better than the second worst position but worse than the best position.
            S(:,1)=[];      S=[S,Vr];
            fit_S(1)=[];    fit_S=[fit_S,fit_Vr];
        elseif fit_Vr<=fit_Vh % Contraction: the reflection performs better than the worst position but worse than the second worst position.
            Vc=Vx+gamma.*(Vr-Vx);
            for d=1:D
                Vc(d)=min(xMax,Vc(d));
                Vc(d)=max(xMin,Vc(d));
            end
            fit_Vc=rotated_rest(Vc',M);%*******************************
            if fit_Vc<fit_Vh % The contruction performs better than the worst position.
                S(:,1)=[];  S=[S,Vc];
                fit_S(1)=[];    fit_S=[fit_S, fit_Vc];
            else % Shrink: the construction performs worse than the worst position.
                for j=1:D
                    S(:,j)=S(:,end)+(S(:,j)-S(:,end)).*sigma; 
                    for d=1:D
                        S(d,j)=max(xMin,S(d,j));
                        S(d,j)=min(xMax,S(d,j));
                    end
                    fit_S(j)=rotated_rest(S(:,j)',M);
                end
            end
        else % Reflection performs worse than the worst position
            Vc=Vx+gamma.*(Vh-Vx);% Sample between the centroid and the worst position.
            for d=1:D
                Vc(d)=max(xMin,Vc(d));
                Vc(d)=min(xMax,Vc(d));
            end
            fit_Vc=rotated_rest(Vc',M);%*********************************
            if fit_Vc<fit_Vh
                S(:,1)=[];  S=[S,Vc];
                fit_S(1)=[];    fit_S=[fit_S, fit_Vc];
            else % Shrink: the new position still performs worse than the worst position.
                for j=1:D
                    S(:,j)=S(:,end)+(S(:,j)-S(:,end)).*sigma; 
                    for d=1:D
                        S(d,j)=max(xMin,S(d,j));
                        S(d,j)=min(xMax,S(d,j));
                    end
                    fit_S(j)=rotated_rest(S(:,j)',M);
                end
            end
        end
        [fit_S,ind]=sort(fit_S,'descend');
        S_tmp=S;
        S=S_tmp(:,ind);
        [yy16(i+1),ind]=min(fit_S);
        Best=S(:,ind);
        if min(fit_S)> des_val
            PDF_cal16 = PDF_cal16 + 1;
        end
   end
   PDF_num_proposed_downhill(iiii) = PDF_cal16;
    if yy16(maxGen - 5) < des_val
        proposed_downhill_16(iiii) = 1;
    end
    Y = 100*ones(1,maxGen+N_add+1);          % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy16;
    ff_end16(iiii,:)=Y;
end
figure
x=1:maxGen+N_add+1;
xlim([0,maxGen+N_add+1]);
semilogy(x,mean(ff_end16,1),'-','linewidth',2,'color',[0.06,1.00,1.00]);
hold on;
xlabel('iterations');
ylabel('fitness');
legend('EVOLER_D','location','SouthWest');
xlim([0,maxGen+N_add])
figure
[~,id_16] = find(proposed_downhill_16 == 1);
PDF_num_proposed_downhill_16 = PDF_num_proposed_downhill(id_16);
[y_pso,x_pso]=hist(PDF_num_proposed_downhill_16,20);
if isempty(PDF_num_proposed_downhill_16)
    PDF_num_proposed_downhill_16 = maxGen;
    p16 = 0;
else
    p16 = length(PDF_num_proposed_downhill_16)/num_all;
end
mean_various=mean(PDF_num_proposed_downhill_16);
Probablity=p16;
h16 = plot(mean_various,Probablity,'>','MarkerSize',12,'lineWidth',1.5,'color',[0.70,0.50,0.70]);
set(h16,'MarkerFaceColor',[0.70,0.50,0.70]);
legend('EVOLER_D')
