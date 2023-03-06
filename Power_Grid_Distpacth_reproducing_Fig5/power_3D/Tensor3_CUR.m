function A_result = Tensor3_CUR(C1_pr,C2_pr,C3_pr,R,c_ll,c_l2,c_l3)
% c = id_all - 1;
% n_all = size(A);
% n1 = n_all(1);
% n2 = n_all(2);
% n3 = n_all(3);
% 
% c_ll = [randperm(n1,c) id_all(1) id_all(1)+1];
% c_l2 = [randperm(n2,c) id_all(1) id_all(1)+1];
% c_l3 = [randperm(n3,c) id_all(1) id_all(1)+1];
% 
% 
% 
% R =tensor(A(c_ll,c_l2,c_l3));
% C1_pr = A(:,c_l2,c_l3);
% C2_pr = A(c_ll,:,c_l3);
% C3_pr = A(c_ll,c_l2,:);


C1 = tenmat(C1_pr,1);
C2 = tenmat(C2_pr,2);
C3 = tenmat(C3_pr,3);


U1 = C1(c_ll,:);
U2 = C2(c_l2,:);
U3 = C3(c_l3,:);


W1 = double(C1*pinv(U1));
W2 = double(C2*pinv(U2));
W3 = double(C3*pinv(U3));


Pro_1 = ttm(R, W1, 1);
Pro_2 = ttm(Pro_1, W2, 2);
Pro_3 = ttm(Pro_2, W3, 3);

A_result = Pro_3;
end