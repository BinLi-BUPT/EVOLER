% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei
close all
clear all

x_all = linspace(-500,500,100);
y_all = linspace(-500,500,100);
%% schewfel: the standard test function
for iii = 1 : length(x_all)
    for jjj = 1 : length(y_all)
        xx = [x_all(iii), y_all(jjj)];
        d = length(xx);
        sum1 = 0;
        for ii = 1:d
            xi = xx(ii);
            sum1 = sum1 + xi*sin(sqrt(abs(xi)));
        end
        y = 418.9829*d - sum1;
        Z(iii,jjj) = y;
    end
end
[m,n] = size(Z);

%% Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples
XXX=[];
for kkkk=1:500
    for jjjj=1:1
        
        % step (i): Structured random sampling on the whole problem space;
        s = 3;    % The sampling length s-O(rlog(r)), r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).  
        Index_c = randi([1,n-10],s,1);
        C1 = Z(:,Index_c);
        Index_r = randi([1,m-10],s,1);
        R1 = Z(Index_r,:);
        Z_samp = zeros(size(Z));
        Z_samp(:,Index_c) = C1;
        Z_samp(Index_r,:) = R1;
        
        % step (ii): Reconstruction of the approximate problem space, via the special form: \hat(Z) = C1 * U1 * R1;
        U0 = Z(Index_r,Index_c);
        [U_u,S_u,V_u]=svd(U0,'econ');
        Uu=U_u';
        s_0 = 3;
        U1 = V_u(:,1:s_0) * pinv(S_u(1:s_0,1:s_0)) * Uu(1:s_0,:);  % Determining the central matrix U1.
        Z_est = C1*U1*R1;                                          % Reconstructing the problem space \hat(Z).
    end
    
    xxx=reshape(Z_est-Z,length(Z)*length(Z),1);                    % The reconstruction error
    XXX = [XXX;xxx];
    kkkk
end

% Figure 1-a£º the reconstruct problem space
figure
mesh(Z_est)
xlabel('x_1');
ylabel('x_2');
zlabel('Reconstructed Function Value');

% The remonstration error between reconstruct problem space and original problem space
figure
mesh(Z-Z_est)
xlabel('x_1');
ylabel('x_2');
zlabel('Reconstructed Error');

% Figure 2-b: the probability density of the residual error of a representation space.
figure
zzz=XXX(find(XXX>-2e-11 & XXX<2e-11));
hist(zzz,60)
xlabel('Residual Error');
ylabel('Histogram of Error');


