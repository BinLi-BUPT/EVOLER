% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;

%%
%-----------initialization parameter setting-------------------------
des_val= 1e-15;                                % If the output result is less than the threshold 'des_val', the global optimum is found  
maxGen = 600;                                  % Maximum number of generations
D= 30;                                         % Dimension
sizepop = 100;                                 % Population size
popMax= 20;                                    % Variable Range [-popmax: popmax]
popMin=-popMax;

s = 2;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).
num_all = 101;                                 % The number of simulation trials
Fitness_all = [];
proposed_downhill_16=zeros(1,num_all);
for iiii = 1 : num_all
    iiii
    %% Step 1: Low-rank Representation Learning based on the dichotomy strategy, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    shift_op = (-2+4*rand(1,D));               % Shifted global optimum, range [-2,2]
    Trail_time = 19;                           % The repeating times of dichotomy strategy based on low-rank sampling
    Discret_leng = 4;                          % Discrete size of the original problem space;      
    max_ite = repmat(popMax,1,D);
    min_ite = repmat(-popMax,1,D);
    tic
    % identification of the global optimum based on low-rank sampling strtegy
    %------ Recon_30D_Initial_shift: Input: Maximum search range, Minimum search range, sampling length, Discrete size
    %---------------------------Output:Global optimum, Maximum search range, Minimum search range
    for ii = 1 : Trail_time
        [global_optimum,pop_div1,pop_div2] = Recon_30D_Initial_shift(max_ite,min_ite,s,Discret_leng,D,shift_op);
        max_ite = pop_div1;
        min_ite = pop_div2;
    end
    N_add = round(Trail_time*((4*Discret_leng*s^4-3*s^5 + Discret_leng^2*s^5)*2+...
        (4*Discret_leng*s^5-3*s^5 + Discret_leng^2*s^4)*2+ (3*Discret_leng*s^5-2*s^5) *2)/sizepop);  % The overall required samples to reconstruct an attention subspace, which is used in the final plot
     
    %% Step 2: Downhill method, which is used to exploit the identified Attention Subspace and finally gain the global optimum.
    %----------------------------------
    % parameter setting
    sigma=0.5; % shrink coefficient
    alpha=1; % reflection coefficient
    beta=2; % expansion coefficient
    gamma=0.5;% contraction coefficient
    xMax=popMax;
    xMin=popMin;
    R=pop_div1-pop_div2;
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
        fit_S(i)=shifted_levy(S(:,i)',shift_op);
    end
    % Sort in descending order of fitness
    [fit_S,ind]=sort(fit_S,'descend');
    S_tmp=S;
    S=S_tmp(:,ind);
    PDF_cal16=N_add;
    yy16(1)=min(fit_S);
    time_1(iiii) = toc;
   for i=1:maxGen
        tic
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
        fit_Vr=shifted_levy(Vr',shift_op);%***********************
        if fit_Vr<=fit_Vl % Expansion:he reflection perform berter than the best position. 
            Ve=Vx+beta.*(Vr-Vx);
            for d=1:D
                Ve(d)=min(xMax,Ve(d));
                Ve(d)=max(xMin,Ve(d));
            end
            fit_Ve=shifted_levy(Ve',shift_op);%******************
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
            fit_Vc=shifted_levy(Vc',shift_op);%*******************************
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
                    fit_S(j)=shifted_levy(S(:,j)',shift_op);
                end
            end
        else % Reflection performs worse than the worst position
            Vc=Vx+gamma.*(Vh-Vx);% Sample between the centroid and the worst position.
            for d=1:D
                Vc(d)=max(xMin,Vc(d));
                Vc(d)=min(xMax,Vc(d));
            end
            fit_Vc=shifted_levy(Vc',shift_op);%*********************************
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
                    fit_S(j)=shifted_levy(S(:,j)',shift_op);
                end
            end
        end
        [fit_S,ind]=sort(fit_S,'descend');
        S_tmp=S;
        S=S_tmp(:,ind);
        [yy16(i),ind]=min(fit_S);
        Best=S(:,ind);
        if min(fit_S)> des_val
            PDF_cal16 = PDF_cal16 + 1;
        end
        time_2(i) = toc;
   end
   PDF_num_proposed_downhill(iiii) = PDF_cal16;
    if yy16(maxGen - 5) < des_val
        proposed_downhill_16(iiii) = 1;
    end
    Y = 100*ones(1,maxGen+N_add);          % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy16;
    Fitness_all(iiii,:)=Y;
    time_evoler_d(iiii) = time_1(iiii)+sum(time_2(1:1200-N_add));
end

% Fit_pro = Fitness_all';
Fit_pro = Fitness_all(2:end,1:900);
T2 = mean(time_2);
save('Downhill_Pro_data_levy','Fit_pro')
save  time_evoler_d1 T2

