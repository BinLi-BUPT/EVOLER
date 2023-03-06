% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clear all;
close all;
%%----------------6-dimensional Power Grid Dispatch--------------%%
% In this code, we transform the original 6-dimensional problem space of into another equivalent 3-dimensional problem,
% i.e., we combining the original 1,2 dimension as one dimension, original 3,4
% dimension as one dimension, original 5, 6 dimension as one dimension
% dimension, which can further reduce the number of samples.

Discret_leng = 11;          
x_end_axis1 = [100:50:500]; % The discrete power range of first generator
x_end_axis2 = [50:10:200];  % The discrete power range of second generator
x_end_axis3 = [80:10:300];  % The discrete power range of third generator
x_end_axis4 = [50:10:150];  % The discrete power range of fourth generator
x_end_axis5 = [50:10:200];  % The discrete power range of fifth generator
x_end_axis6 = [50:10:120];  % The discrete power range of sixth generator
Total = 1260;               % Total power limitations

A_test_ind_en = [];
s = 2;  % the sampling length s-O{rlog(r)}, r is the estimated rank of this problem space (see the online estimation of unknown rank value in Online_rank_estimation.m).                               
t11 = randperm(length(x_end_axis2),s); % randomly sampling index for the first dimension after combining 
t12 = randperm(length(x_end_axis4),s); % randomly sampling index for the second dimension after combining 
t13 = randperm(length(x_end_axis6),s); % randomly sampling index for the third dimension after combining 

% step (i): structured random sampling on the whole 3-dimensional problem space;
% sampling the first small tensor
index_1 = 0;
for iii=1:1:length(x_end_axis1)
    iii;
    for jjj=1:1:length(x_end_axis2)
        index_1 = index_1 + 1;
        index_2 = 0;
        for kkk1 = 1:1:1
        for kkk2 = 1:1:length(t12)
            index_2 = index_2 + 1;
            index_3 = 0;
        for iii1 = 1:1:1
        for iii2 = 1:1:length(t13)
            index_3 = index_3 + 1;
            % function power_allocation_6D: to calculate the cost of 6 generators, the
            % related parameters are provided by ref [12-14] in our manuscript
            % Input: the power of each generator
            % Output: cost.
            A_pro1(index_1,index_2,index_3)=power_allocation_6D([x_end_axis1(iii),x_end_axis2(jjj),x_end_axis3(kkk1),x_end_axis4(t12(kkk2)),x_end_axis5(iii1),x_end_axis6(t13(iii2))]);
        end
        end
        end
        end
    end
end


% sampling the second small tensor
index_1 = 0;
for iii=1:1:1
    iii;
    for jjj=1:1:length(t11)
        index_1 = index_1 + 1;
        index_2 = 0;
        for kkk1 = 1:1:length(x_end_axis3)
        for kkk2 = 1:1:length(x_end_axis4)
            index_2 = index_2 + 1;
            index_3 = 0;
        for iii1 = 1:1:1
        for iii2 = 1:1:length(t13)
            index_3 = index_3 + 1;
            A_pro2(index_1,index_2,index_3)=power_allocation_6D([x_end_axis1(iii),x_end_axis2(t11(jjj)),x_end_axis3(kkk1),x_end_axis4((kkk2)),x_end_axis5(iii1),x_end_axis6(t13(iii2))]);
        end
        end
        end
        end
    end
end

% sampling the third small tensor
index_1 = 0;
for iii=1:1:1
    iii;
    for jjj=1:1:length(t11)
        index_1 = index_1 + 1;
        index_2 = 0;
        for kkk1 = 1:1:1
        for kkk2 = 1:1:length(t12)
            index_2 = index_2 + 1;
            index_3 = 0;
        for iii1 = 1:1:length(x_end_axis5)
        for iii2 = 1:1:length(x_end_axis6)
            index_3 = index_3 + 1;
            A_pro3(index_1,index_2,index_3)=power_allocation_6D([x_end_axis1(iii),x_end_axis2(t11(jjj)),x_end_axis3(kkk1),x_end_axis4(t12(kkk2)),x_end_axis5(iii1),x_end_axis6((iii2))]);
        end
        end
        end
        end
    end
end

% sampling the core tensor
T_ten = [];
index_1 = 0;
for iii=1:1:1
    iii;
    for jjj=1:1:length(t11)
        index_1 = index_1 + 1;
        index_2 = 0;
        for kkk1 = 1:1:1
        for kkk2 = 1:1:length(t12)
            index_2 = index_2 + 1;
            index_3 = 0;
        for iii1 = 1:1:1
        for iii2 = 1:1:length(t13)
            index_3 = index_3 + 1;
            T_ten(index_1,index_2,index_3)=power_allocation_6D([x_end_axis1(iii),x_end_axis2(t11(jjj)),x_end_axis3(kkk1),x_end_axis4(t12(kkk2)),x_end_axis5(iii1),x_end_axis6(t13(iii2))]);
        end
        end
        end
        end
    end
end

C1 = tenmat(tensor(A_pro1),1);  
C2 = tenmat(tensor(A_pro2),2);  
C3 = tenmat(tensor(A_pro3),3);  
Ten_1 = tensor(T_ten);          

% step (ii): reconstruction of the approximate high-dimension problem space
U1 = double(tenmat(Ten_1,1));
U2 = double(tenmat(Ten_1,2));
U3 = double(tenmat(Ten_1,3));
W1 = double(C1*pinv(U1));
W2 = double(C2*pinv(U2));
W3 = double(C3*pinv(U3));
Pro_1 = ttm(Ten_1, W1, 1);
Pro_2 = ttm(Pro_1, W2, 2);
Pro_3 = ttm(Pro_2, W3, 3);
A_result = Pro_3;
A_result2 = A_result;  % the recovered 3-dimensional problem space 

% step (iii): identification of the global optimum in the 3-dimensional representation space, and the determination of attention subspace
index_1 = 1;
x_new_all1 = [];
for ii = 1 : length(x_end_axis1)
    for jj = 1 : length(x_end_axis2)
        index_search1(index_1,:) = [ii,jj]; 
        x_new_all1(index_1) = x_end_axis1(ii) + x_end_axis2(jj);
        index_1 = index_1 + 1;
    end
end

index_1 = 1;
x_new_all2 = [];
for ii = 1 : length(x_end_axis3)
    for jj = 1 : length(x_end_axis4)
        index_search2(index_1,:) = [ii,jj]; 
        x_new_all2(index_1) = x_end_axis3(ii) + x_end_axis4(jj);
        index_1 = index_1 + 1;
    end
end

index_1 = 1;
x_new_all3 = [];
for ii = 1 : length(x_end_axis5)
    for jj = 1 : length(x_end_axis6)
        index_search3(index_1,:) = [ii,jj]; 
        x_new_all3(index_1) = x_end_axis5(ii) + x_end_axis6(jj);
        index_1 = index_1 + 1;
    end
end

A_fina = [];
index_save = [];
for iii=1:1:length(x_new_all1)
    iii;
    for jjj=1:1:length(x_new_all2)
        for kkk1 = 1:1:length(x_new_all3)
            x3 = Total - x_new_all1(iii) - x_new_all2(jjj);
            x3_o = x_new_all3(kkk1);
            if x3_o == x3
                A_fina = [A_fina A_result2(iii,jjj,kkk1)];
                index_save = [index_save;[iii,jjj,kkk1]];
            end
        end
    end
end
[~,index_small] = min(A_fina);
index_save_end = index_save(index_small,:);

% The center of the attention subspace, i.e., the global optimum in the representation space
s_all = [index_search1(index_save_end(1),:), index_search2(index_save_end(2),:), index_search3(index_save_end(3),:)];
P_find = [x_end_axis1(s_all(1)),x_end_axis2(s_all(2)),x_end_axis3(s_all(3)),x_end_axis4(s_all(4)),x_end_axis5(s_all(5)),x_end_axis6(s_all(6))];
save('6D_result','P_find')
