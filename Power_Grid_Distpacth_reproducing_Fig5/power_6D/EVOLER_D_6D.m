% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clear all;
close all;

%-----------initialization parameter setting-------------------------
maxGen = 200;                  % maximum number of generations
sizepop = 20;                % population size
D = 5;                       % dimension
xSum = 1260;
xMax=[500,200,300,150,200];
xMin=[100,50,80,50,50];
DMin=50; DMax=120;
%% Downhill simplex parameter setting
sigma=0.5; % shrink coefficient
alpha=1; % reflection coefficient
beta=2; % expansion coefficient
gamma=0.5;% contraction coefficient
S=zeros(D,D+1);

%(i): load the estimated global optimum in the representation space. Note that, the low-rank represention of original problem space
% has been done by the source code file 'Low_Rank_Represention'.
load('6D_result');     % run the file 'Low_Rank_Represention' to generate this data '6D_result'
x_index_all = P_find;  % the center of the attention subspace, i.e., the global optimum in the representation space
pop = [];
pop2 = [];
N_all=  50;            % the iteration times for performance evaluate
N_add = floor((144*2*2+253*2*2+128*2*2-2*8)/sizepop); % the overall required samples to reconstruct an attention subspace
Y = 1.58e4*ones(N_all,N_add+maxGen);

%(ii): Evolutionary PSO method, which is used to exploit the identified attention subspace and finally gain the global optimum.
% Initialization of population, around the identified attention subspace
for iii = 1 : N_all
    iii
    % Build the initial simplex
    for i=1:(D+1)
        if i==1
            S(:,i)=(x_index_all(1:D)); 
        else
            S(:,i)=S(:,i-1);
            S(i-1,i)=S(i-1,i-1)+4;
        end
        x=S(:,i);
        x(D+1)= xSum-sum(S(:,i));
        fit_S(i)=power_allocation_6D(x);
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
        fit_Vr=power_allocation_6D(xr);
        if fit_Vr<=fit_Vl % Expansion:he reflection perform berter than the best position. 
            Ve=Vx+beta.*(Vx-Vh);
            for d=1:D
                Ve(d)=min(xMax(d),Ve(d));
                Ve(d)=max(xMin(d),Ve(d));
            end
            xe=Ve;
            xe(D+1)=xSum-sum(Ve);
            fit_Ve=power_allocation_6D(xe);
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
            fit_Vc=power_allocation_6D(xc);
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
                    fit_S(j)=power_allocation_6D(x);
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
            fit_Vc=power_allocation_6D(xc);
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
                    fit_S(j)=power_allocation_6D(x);
                end
            end
        end
        [fit_S,ind]=sort(fit_S,'descend');
        S_tmp=S;
        S=S_tmp(:,ind);
        [yy13(i+1),ind]=min(fit_S);
        Best=S(:,ind);
    end
    Y(iii,N_add+1:end) = yy13(1:maxGen);
end
save('downhill_EVOLER_6D','Y')
