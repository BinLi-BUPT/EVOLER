function [Initial_index,p_max,p_min] = Recon_2D_Initial_shift(popmin,popmax,s,Discret_leng,shift_op)

Index_c = randperm(Discret_leng,s); % Randomly column index set
Index_r = randperm(Discret_leng,s); % Randomly row index set
for ii = 1 : length(popmin)
    all_axis(ii,:) = linspace(popmin(ii),popmax(ii),Discret_leng);
end

for iii=1:1:Discret_leng
    for jjj=1:1:length(Index_c)
        C1(iii,jjj)=shifted_weierstrass([all_axis(1, iii),all_axis(2, Index_c(jjj))], shift_op);  % weierstrass: a standard test function
    end
end
for iii=1:1:length(Index_r)
    for jjj=1:1:Discret_leng
        R1(iii,jjj)=shifted_weierstrass([all_axis(1, Index_r(iii)),all_axis(2, (jjj))], shift_op);  % weierstrass: a standard test function
    end
end
U0 = C1(Index_r,:);
U01 = R1(:,Index_c);
%----------------------------------
% step (ii): Reconstruction of the approximate problem space, via the special form: \hat(Z) = C1 * U1 * R1;
% U0 = Z(Index_r,Index_c);
[U_u,S_u,V_u]=svd(U0);
Uu=U_u';
s_0= 2;
U1 = V_u(:,1:s_0) * pinv(S_u(1:s_0,1:s_0)) * Uu(1:s_0,:); % Determining the central matrix U1.
Z_est = C1*U1*R1;                                         % Reconstructing the problem space \hat(Z).

%----------------------------------
% step (iii): Identification of the global optimum in a representation space, and determination of the attention subspace
[m0,n0]=find(Z_est == min(min(Z_est)));                   % The center of the attention subspace, i.e., the global optimum in the representation space
Index_fro = [m0,n0];
for ii = 1 : length(popmin)
    p_max(ii) =all_axis(ii,Index_fro(ii)) +(all_axis(ii,2)-all_axis(ii,1))/2;
    p_min(ii) =all_axis(ii,Index_fro(ii)) -(all_axis(ii,2)-all_axis(ii,1))/2;
    Initial_index(ii) = all_axis(ii,Index_fro(ii));
end
% sum(abs(Initial_index - shift_op))
% Initial_min - shifted_weierstrass(Initial_index,shift_op)

end