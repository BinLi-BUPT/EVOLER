% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
clc;
format long;
close all

%%
%-----------initialization parameter setting-------------------------
des_val= 1e-14;                                % If the output result is less than the threshold 'des_val', the global optimum is found 
maxGen=500;                                    % Maximum number of generations
D=2;                                           % Dimension
N=50;                                          % Population size
popMax= 32;                                    % Variable Range [popmin: popmax]
popMin=-popMax;                                % Variable Range [popmin: popmax]
proposed_downhill_16=zeros(1,maxGen);

Discret_leng = 800;                            % Discrete size of the original problem space; M=N=Discret_leng;
x_axis = linspace(popMin,popMax,Discret_leng); 
y_axis = linspace(popMin,popMax,Discret_leng); 

Discret_leng_now = 1:1:Discret_leng;  
Downsample_factor = 40;
Discret_leng_int = Discret_leng_now(1:Downsample_factor:end); % Down-sampling with the factor 40
Discret_leng_int_num = length(Discret_leng_int);

num_all = 500;  % The number of simulation trials
for iiii = 1 : num_all
    iiii
    % The precomputed original problem space, for the subsequent sampling and reconstruction;
    % This is used to evaluate the reconstruction residual error of the learned low-rank representation.
    M = rot_matrix(D); % rot_matrix:randomly generate one rotated D-dimension matrix: Input: D--dimension; Output: randomly rotated matrix
    for iii=1:1:length(Discret_leng_int)
        for jjj=1:1:length(Discret_leng_int)
            Z(iii,jjj)=ackley_rot([x_axis(Discret_leng_int(iii)),y_axis(Discret_leng_int(jjj))],M); % ackley_rot: a standard test function
        end   
    end
    
   %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples. 
    %----------------------------------
    % step (i): Structured random sampling on the whole problem space,
    % i.e., sampling s rows and columns of the original discrete problem space;
    s = 19;                                         % The sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m). 
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

    N_add = round((s*(Discret_leng_int_num)+s*(Discret_leng_int_num)-s*s)/N); % The overall required samples to reconstruct an attention subspace, which is used in the final plot.  
    
   %% Step 2: Evolutionary downhill method, which is used to exploit the identified Attention Subspace and finally gain the global optimum. 
   %----------------------------------
    % parameter setting
    sigma=0.5; % shrink coefficient
    alpha=1; % reflection coefficient
    beta=2; % expansion coefficient
    gamma=0.5;% contraction coefficient
    x_po1 = x_axis(m0);
    y_po1 = y_axis(n0);
    see2(iiii,:) = [x_po1,y_po1]; 
    n=D+1;
    theat_all = 3*pi*rand(n,1);
    xMax=popMax;
    xMin=popMin;
    % Build the initial simplex 
    for i=1:n
        if i==1
            S(:,i)=[x_axis(m0(1)),y_axis(n0(1))];
        else
            S(:,i)=S(:,i-1);
            S(i-1,i)=S(i-1,i-1)+R;
        end
        fit_S(i)=ackley_rot(S(:,i)',M);
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
        fit_Vr=ackley_rot(Vr',M);%***********************
        if fit_Vr<=fit_Vl % Expansion:he reflection perform berter than the best position. 
            Ve=Vx+beta.*(Vx-Vh);
            for d=1:D
                Ve(d)=min(xMax,Ve(d));
                Ve(d)=max(xMin,Ve(d));
            end
            fit_Ve=ackley_rot(Ve',M);%******************
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
            fit_Vc=ackley_rot(Vc',M);%*******************************
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
                    fit_S(j)=ackley_rot(S(:,j)',M);
                end
            end
        else % Reflection performs worse than the worst position
            Vc=Vx+gamma.*(Vh-Vx);% Sample between the centroid and the worst position.
            for d=1:D
                Vc(d)=max(xMin,Vc(d));
                Vc(d)=min(xMax,Vc(d));
            end
            fit_Vc=ackley_rot(Vc',M);%*********************************
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
                    fit_S(j)=ackley_rot(S(:,j)',M);
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
    Y = 100*ones(1,maxGen+1);          % Assuming the large fitness in the first reconstruction stage
    Y(N_add+1:end) = yy16(1:end-N_add);
    ff_end16(iiii,:)=Y;
end
Fit_pro =  ff_end16;
save('Pro_downhill_data','Fit_pro');
figure()
x=1:maxGen+1;
xlim([0,maxGen]);
semilogy(x,mean(ff_end16,1),'-','linewidth',2,'color',[0.06,1.00,1.00]);
hold on;
xlabel('iterations');
ylabel('fitness');
legend('downhill EVOLER','location','SouthWest');
xlim([0,maxGen])
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
set(h16,'MarkerFaceColor',get(h16,'color'));
legend('downhill EVOLER')
