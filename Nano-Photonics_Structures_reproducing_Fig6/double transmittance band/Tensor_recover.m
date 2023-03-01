function m_all = Tensor_recover(C1_pr,C2_pr,C3_pr,C4_pr,C5_pr,R, C_all)
%% Tensor recovery
C1 = tenmat(C1_pr,1);
C2 = tenmat(C2_pr,2);
C3 = tenmat(C3_pr,3);
C4 = tenmat(C4_pr,4);
C5 = tenmat(C5_pr,5);

c_ll = C_all{1};
c_l2 = C_all{2};
c_l3 = C_all{3};
c_l4 = C_all{4};
c_l5 = C_all{5};

U1 = C1(c_ll,:);
U2 = C2(c_l2,:);
U3 = C3(c_l3,:);
U4 = C4(c_l4,:);
U5 = C5(c_l5,:);
R1 = tucker_als(R, [3,4,4,4,4]); %[4,5,5,5,5]0.96
W1 = double(C1*pinv(U1));
W2 = double(C2*pinv(U2));
W3 = double(C3*pinv(U3));
W4 = double(C4*pinv(U4));
W5 = double(C5*pinv(U5));

Pro_1 = ttm(R1, W1, 1);
Pro_2 = ttm(Pro_1, W2, 2);
Pro_3 = ttm(Pro_2, W3, 3);
Pro_4 = ttm(Pro_3, W4, 4);
Pro_5 = ttm(Pro_4, W5, 5);
A_result = tensor(Pro_5);
[m_all]=find(A_result == collapse(A_result,[1,2,3,4,5],@min));
