% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
clear all;
close all;
clc

%%
%-----------initialization parameter setting-------------------------
xMin=[100,100];
xMax=[600,400];
xSum = 850;                              % total power limitation
maxGen = 400;                            % maximum number of generations
sizepop = 20;                            % population size
D= 2;                
P_i_min = [100,100,50];                  % lower power limitation of each generator
P_i_max = [600,400,200];                 % upper power limitation of each generator
popmax = 500;                            % variable Range [popmin: popmax]
popmin = 50;                             % variable Range [popmin: popmax]


%% Downhill simplex parameter setting
sigma=0.5; % shrink coefficient
alpha=1; % reflection coefficient
beta=2; % expansion coefficient
gamma=0.5;% contraction coefficient
S=zeros(D,D+1);

Discret_leng = 10;                       % discrete size of the original problem space;
x_axis = linspace(popmin,popmax,Discret_leng);
x_axis = x_axis(2:end);
y_axis = linspace(popmin,popmax,Discret_leng);
y_axis = y_axis(2:end);
z_axis = linspace(popmin,popmax,Discret_leng);
z_axis = z_axis(2:end);
Discret_leng2 = length(z_axis);


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

% % Singular values of the problem space (D = 3): Reproducing Figure 4-b
% s11 = tenmat(Z,1);
% [~,sig_1,~] = svd(double(s11),'econ');
% figure
% semilogy(diag(sig_1))
% s11 = tenmat(Z,2);
% [~,sig_1,~] = svd(double(s11),'econ');
% hold on
% semilogy(diag(sig_1))
% s11 = tenmat(Z,3);
% [~,sig_1,~] = svd(double(s11),'econ');
% hold on
% semilogy(diag(sig_1))

Z_ten = tensor(Z);
A_all= [];
N_all = 50;                  % The number of simulation trials
c = 2;                       % the sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).  
N_add = round((length(x_axis)*(c)*3 - 2*(c)*(c)*(c))/sizepop);% the overall required samples to reconstruct an attention subspace, which is used in the final plot.
Y = 8.5e3*ones(N_all,N_add+maxGen);
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
                x3 = xSum - x1 - x2;
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
    z_min2 = xSum - x_min2 - y_min2;
    A_all(sss,:) = [x_min2,y_min2,z_min2]; 
    
    
   %% Step 2: Evolutionary PSO method, which is used to exploit the identified attention subspace and finally gain the global optimum. 
   % Initialization  of population, around the identified attention subspace
   for i=1:(D+1)
        if i==1
            S(:,i)=([x_min2, y_min2]); 
        else
            S(:,i)=S(:,i-1);
            S(i-1,i)=S(i-1,i-1)+2;
        end
        x=S(:,i);
        x(D+1)= xSum-sum(S(:,i));
        fit_S(i)=quad_cf(x);
    end
    % Sort in descending order of fitness
    [fit_S,ind]=sort(fit_S,'descend');
    S_tmp=S;
    S=S_tmp(:,ind);
    PDF_cal13=0;
    yy13(1)=min(fit_S);
    % main loop
    for i=1:maxGen
        Vh=S(:,1);
        fit_Vh=fit_S(1);
        Vl=S(:,end);
        fit_Vl=fit_S(end);   
        Vx=mean(S(:,2:end),2);% centroid
        Vr=Vx+alpha.*(Vx-Vh);% reflection
        for d=1:D
            Vr(d)=min(xMax(d),Vr(d));
            Vr(d)=max(xMin(d),Vr(d));
        end
        xr=Vr;
        xr(D+1)=xSum-sum(Vr);
        fit_Vr=quad_cf(xr);
        if fit_Vr<=fit_Vl % Expansion:he reflection perform berter than the best position. 
            Ve=Vx+beta.*(Vx-Vh);
            for d=1:D
                Ve(d)=min(xMax(d),Ve(d));
                Ve(d)=max(xMin(d),Ve(d));
            end
            xe=Ve;
            xe(D+1)=xSum-sum(Ve);
            fit_Ve=quad_cf(xe);
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
                Vc(d)=min(xMax(d),Vc(d));
                Vc(d)=max(xMin(d),Vc(d));
            end
            xc=Vc;
            xc(D+1)=xSum-sum(Vc);
            fit_Vc=quad_cf(xc);
            if fit_Vc<fit_Vh % The contruction performs better than the worst position.
                S(:,1)=[];  S=[S,Vc];
                fit_S(1)=[];    fit_S=[fit_S, fit_Vc];
            else % Shrink: the construction performs worse than the worst position.
                for j=1:D
                    S(:,j)=S(:,end)+(S(:,j)-S(:,end)).*sigma; 
                    for d=1:D
                        S(d,j)=max(xMin(d),S(d,j));
                        S(d,j)=min(xMax(d),S(d,j));
                    end
                    x=S(:,j);
                    x(D+1)=xSum-sum(S(:,j));
                    fit_S(j)=quad_cf(x);
                end
            end
        else % Reflection performs worse than the worst position
            Vc=Vx+gamma.*(Vh-Vx);% Sample between the centroid and the worst position.
            for d=1:D
                Vc(d)=max(xMin(d),Vc(d));
                Vc(d)=min(xMax(d),Vc(d));
            end
            xc=Vc;
            xc(D+1)=xSum-sum(Vc);
            fit_Vc=quad_cf(xc);
            if fit_Vc<fit_Vh
                S(:,1)=[];  S=[S,Vc];
                fit_S(1)=[];    fit_S=[fit_S, fit_Vc];
            else % Shrink: the new position still performs worse than the worst position.
                for j=1:D
                    S(:,j)=S(:,end)+(S(:,j)-S(:,end)).*sigma; 
                    for d=1:D
                        S(d,j)=max(xMin(d),S(d,j));
                        S(d,j)=min(xMax(d),S(d,j));
                    end
                    x=S(:,j);
                    x(D+1)=xSum-sum(S(:,j));
                    fit_S(j)=quad_cf(x);
                end
            end
        end
        [fit_S,ind]=sort(fit_S,'descend');
        S_tmp=S;
        S=S_tmp(:,ind);
        [yy13(i+1),ind]=min(fit_S);
        Best=S(:,ind);
    end
    Y(sss,N_add+1:end) = yy13(1:maxGen);
end
% figure
% plot([1:size(Y,2)],Y')
% xlabel('Genration')
% ylabel('Fitness')
% Y22 = mean(Y,1);
% hold on
% plot([1:size(Y,2)],Y22,'r')
save('downhill_EVOLER_3D','Y')