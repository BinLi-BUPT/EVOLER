function M=rot_matrix(D)
% rand('state',0)
A=randn(D,D);
[M,~]=qr(A);

end