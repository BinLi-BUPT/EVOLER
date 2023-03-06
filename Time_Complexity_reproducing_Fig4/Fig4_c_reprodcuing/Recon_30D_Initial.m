function [Initial_index,p_max,p_min] = Recon_30D_Initial(popmax,popmin,s,Discret_leng)
      %%%% various division axis for D dimension
      D = length(popmax);
      for ii = 1 : D
          all_axis(ii,:) = linspace(popmin(ii),popmax(ii),Discret_leng);
      end
      %% Step 1: Low-rank Representation Learning, which reconstructs the whole problem space from the very limited samples.
      %----------------------------------
      index_temp = 1;
      Index_sample = [randperm(Discret_leng,s)];               % Randomly index
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
      ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)',all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)',all_axis(13,Index_sample(kk)),...
                  all_axis(14:18,index_temp)', all_axis(19,Index_sample(oo)), all_axis(20:24,index_temp)', all_axis(25,Index_sample(pp)), all_axis(26:30,index_temp)'];
      R_core(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
 end;end;end;end;end;
%%%%-------------- 5-1 -------------------%%%%
%%--------Step 1-----------%%
% step 1_1
for ii = 1 : Discret_leng;for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
    ind_cal = [all_axis(1,ii),all_axis(2, Index_sample(jj)),all_axis(3, Index_sample(kk)),all_axis(4, Index_sample(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
     A51_div_1_1(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end
% step 1_2
for ii = 1 : s;for jj = 1 : Discret_leng;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
      ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2,jj),all_axis(3,Index_sample(kk)),all_axis(4,Index_sample(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
      A51_div_1_2(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end
% step 1_3
for ii = 1 : s;for jj = 1 : s;for kk = 1 : Discret_leng;for oo = 1 : s;for pp = 1 : s
     ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2,Index_sample(jj)),all_axis(3,kk),all_axis(4,Index_sample(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
     A51_div_1_3(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end
% step 1_4
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : Discret_leng;for pp = 1 : s
     ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2,Index_sample(jj)),all_axis(3,Index_sample(kk)),all_axis(4,oo),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
    A51_div_1_4(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end
% step 1_5
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
     ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2,Index_sample(jj)),all_axis(3,Index_sample(kk)),all_axis(4,Index_sample(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
     A51_div_1_R(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end
A_result = Tensor4_CUR(A51_div_1_1,A51_div_1_2,A51_div_1_3,A51_div_1_4,tensor(A51_div_1_R),Index_sample);
A_end_test1 = double(permute(A_result, [4,3,2,1,5]));
A51_div_1_c = reshape(A_end_test1,Discret_leng*Discret_leng*Discret_leng*Discret_leng,s);

for ii = 1 : Discret_leng;for jj = 1 : Discret_leng;for kk = 1 : Discret_leng;for oo = 1 : Discret_leng;for pp = 1 : s
    ind_cal = [all_axis(1,(ii)),all_axis(2,(jj)),all_axis(3,(kk)),all_axis(4,(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
    A51_div_1_end(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal);
end;end;end;end;end;


%%--------Step 2-----------%%

for ii = 1 : s;num1 = 0; for jj = 1 : Discret_leng;for kk = 1 : Discret_leng;for oo = 1 : s;for pp = 1 : s;for tt = 1 : s;for uu = 1 : s;num1 = num1+1;
     ind_cal = [all_axis(1:3,index_temp)', all_axis(4,Index_sample(ii)), all_axis(5,(jj)),all_axis(6,(kk)),all_axis(7,Index_sample(oo)),all_axis(8:12,index_temp)',all_axis(13,Index_sample(pp)),all_axis(14:18,index_temp)',...
           all_axis(19,Index_sample(tt)),all_axis(20:24,index_temp)',all_axis(25,Index_sample(uu)),all_axis(26:30,index_temp)'];
%                    A51_div_2_ca(ii,num1) =  hybrid_func1(ind_cal );
      A51_div_2(ii,jj,kk,oo,pp,tt,uu) = hybrid_func1(ind_cal );
end;end;end;end;end;end;end
A_end_test1 = permute(A51_div_2,[1,7,6,5,4,3,2]);
A51_div_2_r = reshape(A_end_test1,s,Discret_leng*Discret_leng*s*s*s*s);
%  norm(tensor(A51_div_2_ca - A51_div_2_r))
% num = 0;
% for ii = 1 : Discret_leng
%     for jj = 1 : Discret_leng
%         for kk = 1 : Discret_leng
%             for oo = 1 : Discret_leng
%                 num = num + 1;
%                 for pp = 1 : s 
%     ind_cal = [all_axis(1,(ii)),all_axis(2,(jj)),all_axis(3,(kk)),all_axis(4,(oo)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(pp)),all_axis(8:30,index_temp)'];
%     A51_div_1_ca(num,pp) = hybrid_func1(ind_cal );
% end;end;end;end;end

for ii = 1 : s;for jj = 1 : s
        ind_cal = [all_axis(1:3,index_temp)', all_axis(4,Index_sample(ii)),all_axis(5:6,index_temp)',all_axis(7,Index_sample(jj)),all_axis(8:30,index_temp)'];
        A51_div_12_u(ii,jj) = hybrid_func1(ind_cal );
end;end
A51 = A51_div_1_c*pinv(A51_div_12_u)*A51_div_2_r;
L51 = [Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,s,s,s,s];
A2_result_end10 = reshape(A51,L51(4),L51(3),L51(2),L51(1),L51(10),L51(9),L51(8),L51(7),L51(6),L51(5));
A51_tensor = permute(A2_result_end10, [4,3,2,1,10,9,8,7,6,5]);

% for ii = 1: s; for jj = 1; for kk = 1; for oo = 1;for pp = 1
% for ii1 = 1; for jj1 = 1 : s;for kk1 = 1 : s;for oo1 = 1 : s;for pp1 = 1 : s
%             R1(ii,jj1,kk1,oo1,pp1) = A51_tensor(Index_sample(ii),jj,kk,oo,pp,ii1,jj1,kk1,oo1,pp1);
% end;end;end;end;end;end;end;end;end;end
% norm(tensor(R1 - R_core))

%%%%-------------- 5-2 -------------------%%%%
%%--------Step 1-----------%%
% step 2_1
for ii = 1 : s; for jj = 1 : Discret_leng; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
     ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,jj), all_axis(8,Index_sample(kk)),  all_axis(9,Index_sample(oo)),  all_axis(10,Index_sample(pp)), ...
           all_axis(11,Index_sample(tt)), all_axis(12:30,index_temp)'];
      A52_div_1_1(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 2_2
for ii = 1 : s; for jj = 1 : s; for kk = 1 : Discret_leng; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
      ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8,(kk)),  all_axis(9,(oo)),  all_axis(10,Index_sample(pp)), ...
          all_axis(11,Index_sample(tt)), all_axis(12:30,index_temp)'];
       A52_div_1_2(ii,jj,kk,oo,pp,tt) =hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 2_3
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : Discret_leng; for pp = 1 : s; for tt = 1 : s
       ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8,Index_sample(kk)),  all_axis(9,(oo)),  all_axis(10,Index_sample(pp)), ...
           all_axis(11,Index_sample(tt)), all_axis(12:30,index_temp)'];
       A52_div_1_3(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 2_4
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : Discret_leng; for tt = 1 : s
       ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8,Index_sample(kk)),  all_axis(9,Index_sample(oo)),  all_axis(10,(pp)), ...
          all_axis(11,Index_sample(tt)), all_axis(12:30,index_temp)'];
        A52_div_1_4(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 2_5
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8,Index_sample(kk)),  all_axis(9,Index_sample(oo)),  all_axis(10,Index_sample(pp)), ...
           all_axis(11,Index_sample(tt)),all_axis(12:30,index_temp)'];
         A52_div_1_R(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
A_result = Tensor4_CUR_52(A52_div_1_1,A52_div_1_2,A52_div_1_3,A52_div_1_4,tensor(A52_div_1_R),Index_sample);
A_end_test1 = double(permute(A_result, [5,4,3,2,1,6]));
A52_div_1_c = reshape(A_end_test1,s*Discret_leng*Discret_leng*Discret_leng*Discret_leng,s);
%%--------Step 2-----------%%
for ii = 1 : s; for jj = 1 : Discret_leng;for kk = 1 : Discret_leng;for oo = 1 : s;for pp = 1 : s;for tt = 1 : s
       ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:10,index_temp)', all_axis(11,(jj)),all_axis(12,(kk)),all_axis(13,Index_sample(oo)),all_axis(14:18,index_temp)',all_axis(19,Index_sample(pp)),...
             all_axis(20:24,index_temp)',all_axis(25,Index_sample(tt)),all_axis(26:30,index_temp)'];
        A52_div_2(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
A_end_test1 = permute(A52_div_2,[1,6,5,4,3,2]);
A52_div_2_r = reshape(A_end_test1,s,Discret_leng*Discret_leng*s*s*s);
% norm(tensor(A51_div_2_ca - A51_div_2_r))
for ii = 1 : s;for jj = 1 : s
        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:10,index_temp)',all_axis(11,Index_sample(jj)), all_axis(12:30,index_temp)'];
        A52_div_12_u(ii,jj) = hybrid_func1(ind_cal );
end;end

A52 = A52_div_1_c*pinv(A52_div_12_u)*A52_div_2_r;
% %% test the accuracy
L52 = [s,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,s,s,s];
A2_result_end10 = reshape(A52,L52(5),L52(4),L52(3),L52(2),L52(1),L52(10),L52(9),L52(8),L52(7),L52(6));
A52_tensor = permute(A2_result_end10, [5,4,3,2,1,10,9,8,7,6]);

% for ii = 1: s; for jj = 1:s; for kk = 1; for oo = 1;for pp = 1
%                     for ii1 = 1; for jj1 = 1;for kk1 = 1 : s;for oo1 = 1 : s;for pp1 = 1 : s
%                                         R2(ii,jj,kk1,oo1,pp1) = A52_tensor((ii),Index_sample(jj),kk,oo,pp,ii1,jj1,kk1,oo1,pp1);
%                     end;end;end;end;end;end;end;end;end;end
% norm(tensor(R2 - R_core))

%%%%-------------- 5-3 -------------------%%%%
%%--------Step 1-----------%%
% step 3_1
for ii = 1 : s; for jj = 1 : s; for kk = 1 : Discret_leng; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,kk), all_axis(14,Index_sample(oo)), all_axis(15,Index_sample(pp)),...
                            all_axis(16,Index_sample(tt)),all_axis(17:30,index_temp)'];
                        A53_div_1_1(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_2
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : Discret_leng; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)), all_axis(14,(oo)), all_axis(15,Index_sample(pp)),...
                            all_axis(16,Index_sample(tt)),all_axis(17:30,index_temp)'];
                        A53_div_1_2(ii,jj,kk,oo,pp,tt) =  hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_3
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : Discret_leng; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)), all_axis(14,Index_sample(oo)), all_axis(15,(pp)),...
                            all_axis(16,Index_sample(tt)),all_axis(17:30,index_temp)'];
                        A53_div_1_3(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_4
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)), all_axis(14,Index_sample(oo)), all_axis(15,Index_sample(pp)),...
                            all_axis(16,Index_sample(tt)),all_axis(17:30,index_temp)'];
                        A53_div_1_R(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end

A_result = Tensor4_CUR_53(A53_div_1_1,A53_div_1_2,A53_div_1_3,tensor(A53_div_1_R),Index_sample);
A_end_test1 = double(permute(A_result, [5,4,3,2,1,6]));
A53_div_1_c = reshape(A_end_test1,s*Discret_leng*Discret_leng*Discret_leng*s,s);

%%--------Step 2-----------%%
% step 3_2_1
for ii = 1 : s; for jj = 1 : Discret_leng; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:15,index_temp)', all_axis(16,(jj)),all_axis(17,Index_sample(kk)), all_axis(18,Index_sample(oo)), all_axis(19,Index_sample(pp)),...
                            all_axis(20:24,index_temp)',all_axis(25,Index_sample(tt)),all_axis(26:30,index_temp)'];
                        A53_div_2_1(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_2_2
for ii = 1 : s; for jj = 1 : s; for kk = 1 : Discret_leng; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:15,index_temp)', all_axis(16,Index_sample(jj)),all_axis(17,(kk)), all_axis(18,Index_sample(oo)), all_axis(19,Index_sample(pp)),...
                            all_axis(20:24,index_temp)',all_axis(25,Index_sample(tt)),all_axis(26:30,index_temp)'];
                        A53_div_2_2(ii,jj,kk,oo,pp,tt) =  hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_2_3
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : Discret_leng; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:15,index_temp)', all_axis(16,Index_sample(jj)),all_axis(17,Index_sample(kk)), all_axis(18,(oo)), all_axis(19,Index_sample(pp)),...
                            all_axis(20:24,index_temp)',all_axis(25,Index_sample(tt)),all_axis(26:30,index_temp)'];
                        A53_div_2_3(ii,jj,kk,oo,pp,tt) =  hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 3_2_4
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:15,index_temp)', all_axis(16,Index_sample(jj)),all_axis(17,Index_sample(kk)), all_axis(18,Index_sample(oo)), all_axis(19,Index_sample(pp)),...
                            all_axis(20:24,index_temp)',all_axis(25,Index_sample(tt)),all_axis(26:30,index_temp)'];
                        A53_div_2_R(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end

A_result = Tensor4_CUR_53_2(A53_div_2_1,A53_div_2_2,A53_div_2_3,tensor(A53_div_2_R),Index_sample);
A_end_test1 = double(permute(A_result, [1 6 5 4 3 2]));
A53_div_1_r = reshape(A_end_test1,s,s*Discret_leng*Discret_leng*Discret_leng*s);
for ii = 1 : s;for jj = 1 : s
        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:15,index_temp)',all_axis(16,Index_sample(jj)), all_axis(17:30,index_temp)'];
        A53_div_2_u(ii,jj) = hybrid_func1(ind_cal );
end;end
A53 = A53_div_1_c*pinv(A53_div_2_u)*A53_div_1_r;
L53 = [s,s,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,s,s];
A2_result_end10 = reshape(A53,L53(5),L53(4),L53(3),L53(2),L53(1),L53(10),L53(9),L53(8),L53(7),L53(6));
A53_tensor = permute(A2_result_end10, [5,4,3,2,1,10,9,8,7,6]);


% for ii = 1: s; for jj = 1:s; for kk = 1:s; for oo = 1;for pp = 1
%                     for ii1 = 1; for jj1 = 1;for kk1 = 1;for oo1 = 1 : s;for pp1 = 1 : s
%                                         R3(ii,jj,kk,oo1,pp1) = A53_tensor((ii),(jj),Index_sample(kk),oo,pp,ii1,jj1,kk1,oo1,pp1);
%                     end;end;end;end;end;end;end;end;end;end
% norm(tensor(R3 - R_core))

%%%%-------------- 5-4 -------------------%%%%
%%--------Step 1-----------%%
num = 0;
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : Discret_leng; for pp = 1 : Discret_leng; num = num + 1; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)),all_axis(8:12,index_temp)', all_axis(13,Index_sample(kk)),...
                            all_axis(14:18,index_temp)',all_axis(19,oo), all_axis(20,pp), all_axis(21:24,index_temp)',all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                        A54_div_1_c(num,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end

%%--------Step 2-----------%%
% step 4_1
for ii = 1 : s; for jj = 1 : Discret_leng; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:20,index_temp)', all_axis(21,(jj)),all_axis(22,Index_sample(kk)),...
                            all_axis(23,Index_sample(oo)), all_axis(24,Index_sample(pp)),all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                        A54_div_2_1(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 4_2
for ii = 1 : s; for jj = 1 : s; for kk = 1 : Discret_leng; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:20,index_temp)', all_axis(21,Index_sample(jj)),all_axis(22,(kk)),...
                            all_axis(23,Index_sample(oo)), all_axis(24,Index_sample(pp)),all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                        A54_div_2_2(ii,jj,kk,oo,pp,tt) =hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 4_3
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : Discret_leng; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:20,index_temp)', all_axis(21,Index_sample(jj)),all_axis(22,Index_sample(kk)),...
                            all_axis(23,(oo)), all_axis(24,Index_sample(pp)),all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                        A54_div_2_3(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 4_4
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : Discret_leng; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:20,index_temp)', all_axis(21,Index_sample(jj)),all_axis(22,Index_sample(kk)),...
                            all_axis(23,Index_sample(oo)), all_axis(24,(pp)),all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                        A54_div_2_4(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
% step 4_5
for ii = 1 : s; for jj = 1 : s; for kk = 1 : s; for oo = 1 : s; for pp = 1 : s; for tt = 1 : s
                        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:20,index_temp)', all_axis(21,Index_sample(jj)),all_axis(22,Index_sample(kk)),...
                            all_axis(23,Index_sample(oo)), all_axis(24,Index_sample(pp)),all_axis(25,Index_sample(tt)), all_axis(26:30,index_temp)'];
                         A54_div_2_R(ii,jj,kk,oo,pp,tt) = hybrid_func1(ind_cal );
end;end;end;end;end;end
A_result = Tensor4_CUR_52(A54_div_2_1,A54_div_2_2,A54_div_2_3,A54_div_2_4,tensor(A54_div_2_R),Index_sample);
A_end_test1 = double(permute(A_result, [1 6 5 4 3 2]));
A54_div_1_r = reshape(A_end_test1,s,Discret_leng*Discret_leng*Discret_leng*Discret_leng*s);
%%% combine
for ii = 1 : s; for jj = 1 : s
        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:24,index_temp)',all_axis(25,Index_sample(jj)), all_axis(26:30,index_temp)'];
        A54_div_2_u(ii,jj) = hybrid_func1(ind_cal );
end;end
A54 = A54_div_1_c*pinv(A54_div_2_u)*A54_div_1_r;

% %% test the accuracy
L54 = [s,s,s,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,s];
A2_result_end10 = reshape(A54,L54(5),L54(4),L54(3),L54(2),L54(1),L54(10),L54(9),L54(8),L54(7),L54(6));
A54_tensor = permute(A2_result_end10, [5,4,3,2,1,10,9,8,7,6]);

% for ii = 1: s; for jj = 1:s; for kk = 1:s; for oo = 1:s;for pp = 1
%                     for ii1 = 1; for jj1 = 1;for kk1 = 1;for oo1 = 1;for pp1 = 1 : s
%                                         R4(ii,jj,kk,oo,pp1) = A54_tensor((ii),(jj),(kk),Index_sample(oo),pp,ii1,jj1,kk1,oo1,pp1);
%                     end;end;end;end;end;end;end;end;end;end
% norm(tensor(R4 - R_core))


%%%%-------------- 5-5 -------------------%%%%
%%--------Step 1-----------%%
num = 0;
for ii = 1 : s; for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : Discret_leng;for tt = 1 : Discret_leng; num = num + 1; for uu = 1 : s
      ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:6,index_temp)', all_axis(7,Index_sample(jj)), all_axis(8:12,index_temp)',all_axis(13,Index_sample(kk)), all_axis(14:18,index_temp)',...
            all_axis(19,Index_sample(oo)),all_axis(20:24,index_temp)',all_axis(25,(pp)),all_axis(26,(tt)),all_axis(27,Index_sample(uu)),all_axis(28:30,index_temp)'];
      A55_div_c(num,uu) = hybrid_func1(ind_cal );
                            %     A51_div_2_ca(ii,num1) =  schewfel(ind_cal);
end;end;end;end;end;end;end

%%--------Step 2-----------%%
%--5-5-1
for ii = 1 : s;for jj = 1 : Discret_leng;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
                    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:26,index_temp)',all_axis(27,(jj)),all_axis(28,Index_sample(kk)),all_axis(29,Index_sample(oo)),all_axis(30,Index_sample(pp))];
                    A55_div_1_1(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal );
end;end;end;end;end
%--5-5-2
for ii = 1 : s;for jj = 1 : s;for kk = 1 : Discret_leng;for oo = 1 : s;for pp = 1 : s
                    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:26,index_temp)',all_axis(27,Index_sample(jj)),all_axis(28,(kk)),all_axis(29,Index_sample(oo)),all_axis(30,Index_sample(pp))];
                    A55_div_1_2(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal );
end;end;end;end;end
%--5-5-3
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : Discret_leng;for pp = 1 : s
                    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:26,index_temp)',all_axis(27,Index_sample(jj)),all_axis(28,Index_sample(kk)),all_axis(29,(oo)),all_axis(30,Index_sample(pp))];
                    A55_div_1_3(ii,jj,kk,oo,pp) = hybrid_func1(ind_cal );
end;end;end;end;end
%--5-5-4
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : Discret_leng
                    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:26,index_temp)',all_axis(27,Index_sample(jj)),all_axis(28,Index_sample(kk)),all_axis(29,Index_sample(oo)),all_axis(30,(pp))];
                    A55_div_1_4(ii,jj,kk,oo,pp) =hybrid_func1(ind_cal );
end;end;end;end;end
%--5-5-5
for ii = 1 : s;for jj = 1 : s;for kk = 1 : s;for oo = 1 : s;for pp = 1 : s
                    ind_cal = [all_axis(1,Index_sample(ii)),all_axis(2:26,index_temp)',all_axis(27,Index_sample(jj)),all_axis(28,Index_sample(kk)),all_axis(29,Index_sample(oo)),all_axis(30,Index_sample(pp))];
                    A55_div_1_R(ii,jj,kk,oo,pp) =hybrid_func1(ind_cal );
end;end;end;end;end

A_result = Tensor4_CUR_52(A55_div_1_1,A55_div_1_2,A55_div_1_3,A55_div_1_4,tensor(A55_div_1_R ),Index_sample);
A_end_test1 = permute(A_result,[1,5,4,3,2]);
A55_div_2_r = reshape(double(A_end_test1),s,Discret_leng*Discret_leng*Discret_leng*Discret_leng);
%%% combine
for ii = 1 : s; for jj = 1 : s
        ind_cal = [all_axis(1,Index_sample(ii)), all_axis(2:26,index_temp)',all_axis(27,Index_sample(jj)), all_axis(28:30,index_temp)'];
        A55_div_2_u(ii,jj) = hybrid_func1(ind_cal );
end;end
A55 = A55_div_c*pinv(A55_div_2_u)*A55_div_2_r;

% %% test the accuracy
L55 = [s,s,s,s,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng,Discret_leng];
A2_result_end10 = reshape(A55,L55(6),L55(5),L55(4),L55(3),L55(2),L55(1),L55(10),L55(9),L55(8),L55(7));
A55_tensor = permute(A2_result_end10, [6,5,4,3,2,1,10,9,8,7]);

% for ii = 1: s; for jj = 1:s; for kk = 1:s; for oo = 1:s;for pp = 1:s
% for ii1 = 1; for jj1 = 1;for kk1 = 1;for oo1 = 1;for pp1 = 1
%      R5(ii,jj,kk,oo,pp) = A55_tensor((ii),(jj),(kk),(oo),Index_sample(pp),ii1,jj1,kk1,oo1,pp1);
%  end;end;end;end;end;end;end;end;end;end
%  norm(tensor(R5 - R_core))
[m_all_1]=find(tensor(A51_tensor) == collapse(tensor(A51_tensor),[1,2,3,4,5,6,7,8,9,10],@min));
[m_all_2]=find(tensor(A52_tensor) == collapse(tensor(A52_tensor),[1,2,3,4,5,6,7,8,9,10],@min));
[m_all_3]=find(tensor(A53_tensor) == collapse(tensor(A53_tensor),[1,2,3,4,5,6,7,8,9,10],@min));
[m_all_4]=find(tensor(A54_tensor) == collapse(tensor(A54_tensor),[1,2,3,4,5,6,7,8,9,10],@min));
[m_all_5]=find(tensor(A55_tensor) == collapse(tensor(A55_tensor),[1,2,3,4,5,6,7,8,9,10],@min));
Index_fro = [m_all_1(1,1:6),m_all_2(1,2:7),m_all_3(1,3:8),m_all_4(1,4:9),m_all_5(1,5:10)];
Initial_min = double(Result_cal(A51_tensor,A52_tensor,A53_tensor,A54_tensor,A55_tensor,R_core,s,Index_fro));
for ii = 1 : D
    p_max(ii) =all_axis(ii,Index_fro(ii)) +(all_axis(ii,2)-all_axis(ii,1))/2;
    p_min(ii) =all_axis(ii,Index_fro(ii)) -(all_axis(ii,2)-all_axis(ii,1))/2;
    Initial_index(ii) = all_axis(ii,Index_fro(ii));
end
% sum(abs(Initial_index - shift_op))
% Initial_min - hybrid_func1(Initial_index)

end