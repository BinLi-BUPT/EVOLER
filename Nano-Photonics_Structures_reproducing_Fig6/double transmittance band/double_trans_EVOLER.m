% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clc
close all
clear all

% To minimize the complex interaction between MATLAB and LUMERICAL, we provide the data of 5-dimensional tensor--'double_trans_results', 
% which directly gives the fitness on a discrete problem space. This was 
% obtained by running the LUMERICAL software in advance, and serves as the
% simple fitness function in the first stage of structured random sampling.
% One this basis, EVOLER can directly determine the sampling indexes and evaluate the fitness, without frequently switching to the LUMERICAL running in the sampling process. 
A_end = load('double_trans_results');             % D=5 layer paramters; a1 \in [15:15:150]nm, a2 \in [15:15:150]nm, a3 \in [15:15:150]nm, a4 \in [15:15:150]nm, a5 \in [15:15:150]nm;
A_ten_end = tensor(A_end.data_rec_1);

%% Low-rank Representation Learning, which reconstructs the whole high-dimensional problem space from the very limited samples
%----------------------------------------------------
% (i) Determine the structured sampling indexes of the whole
% high-dimensional problem space; the detail informaiton of an iterative determination of random
% sampling index can be found in SI Section I.
data_rec_ten = tensor(A_ten_end);
C_sapm = [1,10];       % The initial sampling index, 
A = data_rec_ten;
A_end_end = [];
tao = 0.91;            % Threshold (>0.9) of accumulated energy of the first r singular values, to estimate the rank of problem space and configure a sampling length.

% 1) Determine the sampling index information at first dimension, according
% to the provided small tensor 'A1', which is composed by some elements of
% the original problem space, and it is obtained by running the LUMERICAL
% software in advance.
A1 = A(:,C_sapm,C_sapm,C_sapm,C_sapm);  
R_tem1 = tenmat(A1,1);

% Determine the sampling length of the first dimension according to the singular value of A1
% rat_cal: Input--(small tensor, singular threshold);
%          Output--sampling length
rat_sig1 = rat_cal(R_tem1,tao);  

% Determine the sampling index of the first dimension
% inde_samd: Input--(small tensor, initial sampling index, sampling length);
%            Output--final sampling index of first dimension
c1_temp = inde_samd(R_tem1,C_sapm,rat_sig1+1); 
C_all{1} = c1_temp;

% 2) Determine the sampling index information at second dimension, according
% to the provided small tensor 'A2', which is composed by some elements of
% the original problem space, and it is obtained by running the LUMERICAL
% software in advance.
A2 = A(c1_temp,:,C_sapm,C_sapm,C_sapm);  
R_tem2 = tenmat(A2,2); 

% Determine the sampling length of the second dimension according to the singular value of A2
% rat_cal: Input--(small tensor, singular threshold);
%          Output--sampling length
rat_sig2 = rat_cal(R_tem2,tao);

% Determine the sampling index of the second dimension
% inde_samd: Input--(small tensor, initial sampling index, sampling length);
%            Output--final sampling index of second dimension
c2_temp = inde_samd(R_tem2,C_sapm,rat_sig2);
C_all{2} = c2_temp;


% 3) Determine the sampling index information at third dimension, according
% to the provided small tensor 'A3', which is composed by some elements of
% the original problem space, and it is obtained by running the LUMERICAL
% software in advance.
A3 = A(c1_temp,c2_temp,:,C_sapm,C_sapm);
R_tem3 = tenmat(A3,3); %

% Determine the sampling length of the third dimension according to the singular value of A3
% rat_cal: Input--(small tensor, singular threshold);
%          Output--sampling length
rat_sig3 = rat_cal(R_tem3,tao);

% Determine the sampling index of the third dimension
% inde_samd: Input--(small tensor, initial sampling index, sampling length);
%            Output--final sampling index of third dimension
c3_temp = inde_samd(R_tem3,C_sapm,rat_sig3);
C_all{3} = c3_temp;


% 4) Determine the sampling index information at fourth dimension, according
% to the provided small tensor 'A4', which is composed by some elements of
% the original problem space, and it is obtained by running the LUMERICAL
% software in advance.
A4 = A(c1_temp,c2_temp,c3_temp,:,C_sapm); 
R_tem4 = tenmat(A4,4); 

% Determine the sampling length of the fourth dimension according to the singular value of A4
% rat_cal: Input--(small tensor, singular threshold);
%          Output--sampling length
rat_sig4 = rat_cal(R_tem4,tao);

% Determine the sampling index of the fourth dimension
% inde_samd: Input--(small tensor, initial sampling index, sampling length);
%            Output--final sampling index of fourth dimension
c4_temp = inde_samd(R_tem4,C_sapm,rat_sig4);
C_all{4} = c4_temp;

% 5) Determine the sampling index information at fifth dimension, according
% to the provided small tensor 'A5', which is composed by some elements of
% the original problem space, and it is obtained by running the LUMERICAL
% software in advance.
A5 = A(c1_temp,c2_temp,c3_temp,c4_temp,:);
R_tem5 = tenmat(A5,5); 

% Determine the sampling length of the fifth dimension according to the singular value of A5
% rat_cal: Input--(small tensor, singular threshold);
%          Output--sampling length
rat_sig5 = rat_cal(R_tem5,tao);

% Determine the sampling index of the fifth dimension
% inde_samd: Input--(small tensor, initial sampling index, sampling length);
%            Output--final sampling index of fifth dimension
c5_temp = inde_samd(R_tem5,C_sapm,rat_sig5);
C_all{5} = c5_temp;

% (ii) Reconstruction of the approximate problem space with limited samples
C1_pr = A(:,c2_temp,c3_temp,c4_temp,c5_temp);          % The first sampling small tensor
C2_pr = A(c1_temp,:,c3_temp,c4_temp,c5_temp);          % The second sampling small tensor
C3_pr = A(c1_temp,c2_temp,:,c4_temp,c5_temp);          % The third sampling small tensor
C4_pr = A(c1_temp,c2_temp,c3_temp,:,c5_temp);          % The fourth sampling small tensor
C5_pr = A(c1_temp,c2_temp,c3_temp,c4_temp,:);          % The fifth sampling small tensor
R =tensor(A(c1_temp,c2_temp,c3_temp,c4_temp,c5_temp)); % The sampling core tensor

% Identification of the global optimum in the representation space, the determination of attention subspace 
% Tensor_recover: Input--(sampling small tensor, core tensor, sampling index set);
%                 Output--the center of attention subspace
m_all = Tensor_recover(C1_pr,C2_pr,C3_pr,C4_pr,C5_pr,R, C_all);

T_11 = [15:15:150];
index_2_est = flip(T_11(m_all));                       % the center of attention subspace used for LUMERICAL simulations
num_1 = prod(size(C1_pr));
num_2 = prod(size(C2_pr));
num_3 = prod(size(C3_pr));
num_4 = prod(size(C4_pr));
num_5 = prod(size(C5_pr));
num_52 = prod(size(R));
N_add = num_1+num_2+num_3+num_4+num_5-4*num_52;        % The overall required samples to reconstruct an attention subspace. 
save('index_2_est','index_2_est')                      % Save the estimated global optimum.