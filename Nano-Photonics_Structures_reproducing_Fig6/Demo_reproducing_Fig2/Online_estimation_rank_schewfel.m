% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

% Rank estimation demo for schewfel function

close all
clear all

popmax = 500;                                                % Parameter range of each dimension
popmin = -500;                                               % Parameter range of each dimension


Discret_leng = 20;                                           % Discrete size of each dimension of problems space
x_axis = linspace(popmin,popmax,Discret_leng);
y_axis = linspace(popmin,popmax,Discret_leng);
for iii=1:1:Discret_leng
    for jjj=1:1:Discret_leng
        Z(iii,jjj)=schewfel([x_axis(iii),y_axis(jjj)]);      % Evalutation of test function, schewfel function
    end
end
[m,n] = size(Z);
figure
[~,SS,~]=svd(Z);
semilogy(diag(SS))
xlabel('Index')
ylabel('Singular Values')
hold on;

% mesh(X,Y,Z)
tao = 0.999;   % Threshold of accumulated energy of the first r singular values, to estimate the rank of problem space and configure a sampling length.
XXX=[];
for kkk = 1 : 100
    for iiii = 3:1:20
        s = iiii*1;
        Index_c = randi([1,n-10],s,1);
        C1 = Z(:,Index_c);
        Index_r = randi([1,m-10],s,1);
        R1 = Z(Index_r,:);
        U0 = Z(Index_r,Index_c);
        [U_u,S_u,V_u]=svd(U0,'econ');
        Uu=U_u';
        s_0 = s;
        U1 = V_u(:,1:s_0) * pinv(S_u(1:s_0,1:s_0)) * Uu(1:s_0,:);
        Z_est = C1*(pinv(U0)+U1)*R1/2;
        [~,S_est,~]=svd(Z_est);
        
        singular_valur_est = diag(S_est);
        singular_valur_est = singular_valur_est(1:s-1);
        p1 = sum(singular_valur_est(1:end-1))/sum(singular_valur_est);
        singValue = singular_valur_est;
        %     figure
        %     semilogy(singValue)
        %     xlabel('Index')
        %     ylabel('Singular Values')
        if p1 > tao
            rank_est = s-1;
            break
        else
        end
    end
    rank_est_all(kkk) = rank_est-1;
end
rank_est_end = ceil(mean(rank_est_all))

