function A_result = Tensor4_CUR_53(C1_pr,C2_pr,C3_pr,R,c_ll)
C1 = tenmat(C1_pr,3);
C2 = tenmat(C2_pr,4);
C3 = tenmat(C3_pr,5);
U1 = C1(c_ll,:);
U2 = C2(c_ll,:);
U3 = C3(c_ll,:);
W1 = double(C1*pinv(U1));
W2 = double(C2*pinv(U2));
W3 = double(C3*pinv(U3));
Pro_1 = ttm(R, W1, 3);
Pro_2 = ttm(Pro_1, W2, 4);
Pro_3 = ttm(Pro_2, W3, 5);
A_result = Pro_3;
end