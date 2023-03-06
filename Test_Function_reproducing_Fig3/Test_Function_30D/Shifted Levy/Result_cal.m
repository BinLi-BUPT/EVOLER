function A_result1 = Result_cal(A51_tensor,A52_tensor,A53_tensor,A54_tensor,A55_tensor,R_core,s,test_ind)

a1 =  (A51_tensor(test_ind(1),test_ind(2),test_ind(3),test_ind(4),test_ind(5),test_ind(6),:,:,:,:));
C1_pr = (reshape(a1,1,s,s,s,s));
a2 = A52_tensor(:,test_ind(7),test_ind(8),test_ind(9),test_ind(10),test_ind(11),test_ind(12),:,:,:);
C2_pr = tensor(reshape(a2,s,1,s,s,s));
a3 = A53_tensor(:,:,test_ind(13),test_ind(14),test_ind(15),test_ind(16),test_ind(17),test_ind(18),:,:);
C3_pr = tensor(reshape(a3,s,s,1,s,s));
a4 = A54_tensor(:,:,:,test_ind(19),test_ind(20),test_ind(21),test_ind(22),test_ind(23),test_ind(24),:);
C4_pr = tensor(reshape(a4,s,s,s,1,s));
a5 = A55_tensor(:,:,:,:,test_ind(25),test_ind(26),test_ind(27),test_ind(28),test_ind(29),test_ind(30));
C5_pr = tensor(reshape(a5,s,s,s,s,1));

C1 = tenmat(C1_pr,1);
C2 = tenmat(C2_pr,2);
C3 = tenmat(C3_pr,3);
C4 = tenmat(C4_pr,4);
C5 = reshape(double(C5_pr),1,prod(size(C5_pr)));
U1 = double(tenmat(R_core,1));
U2 = double(tenmat(R_core,2));
U3 = double(tenmat(R_core,3));
U4 = double(tenmat(R_core,4));
U5 = double(tenmat(R_core,5));
W1 = double(C1*pinv(U1));
W2 = double(C2*pinv(U2));
W3 = double(C3*pinv(U3));
W4 = double(C4*pinv(U4));
W5 = double(C5*pinv(U5));
Pro_1 = ttm(tensor(R_core), W1, 1);
Pro_2 = ttm(Pro_1, W2, 2);
Pro_3 = ttm(Pro_2, W3, 3);
Pro_4 = ttm(Pro_3, W4, 4);
Pro_5 = ttm(Pro_4, W5, 5);
A_result1 = Pro_5;
end
