% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clear all
clc
close all;

%% Running several files to generate the following data 
data_1 = load('CPSO_H_6D'); % run the MATLAB file 'CPSO_H_6D.m' to generate this data
data_2 = load('CPSO_S_6D'); % run the MATLAB file 'CPSO_H_6D.m' to generate this data
data_3 = load('DPSO_6D');   % run the MATLAB file 'DPSO_6D.m' to generate this data
data_4 = load('MCPSO_6D');   % run the MATLAB file 'MPSO_6D.m' to generate this data
data_5 = load('PSO_6D');    % run the MATLAB file 'PSO_6D.m' to generate this data
data_6 = load('CLPSO_6D');  % run the MATLAB file 'CLPSO_6D.m' to generate this data
data_7 = load('EVOLER_6D');         % first run the file 'Low_Rank_Represention', and then run the MATLAB file 'EVOLER_6D.m' to generate this data
data_8 = load('NPSO_6D');    % run the MATLAB file 'NPSO_6D.m' to generate this data
data_9 = load('GC_PSO_6D');    % run the MATLAB file 'GC_PSO_6D.m' to generate this data
data_10 = load('DNSPSO_6D');    % run the MATLAB file 'DNS_PSO_6D.m' to generate this data
data_11 = load('MPCPSO_6D');   % run the MATLAB file 'MPC_PSO_6D.m' to generate this data
data_12 = load('LPSO_6D');     %run the MATLAB file 'LPSO_6D.m' to generate this data
data_14 = load('LPSO_EVOLER_6D');       %run the MATLAB file 'EVOLER_L_6D.m' to generate this data
data_15 = load('downhill_EVOLER_6D');    %run the MATLAB file 'EVOLER_D_6D.m' to generate this data


data_CPSO_H = data_1.COSTS;
data_CPSO_S = data_2.COSTS;
data_DPSO = data_3.COSTS;
data_MPSO = data_4.COSTS;
data_PSO = data_5.COSTS;
data_CLPSO = data_6.COSTS;
data_proposed = data_7.Y;
data_NPSO = data_8.Y;
data_GC_PSOd = data_9.Y;
data_DNSPSO = data_10.Y;
data_MPCPSO = data_11.Y;
data_LPSO=data_12.COSTS;
data_LPSO_EVOLER=data_14.Y;
data_downhill_EVOLER=data_15.Y;

limit_de = 1e-3;
% the decision threshold Tao_1, if the PSO searching result is less than
% this threshold, the minimum has been found.
Tao_1 = (1.52915e+04)*(1+0.001);

% To search the number of times that different comparing method can find the minimum value in 50 iterations 
% and what they are respectively
%% 1
diff_CPSO_H = diff(data_CPSO_H',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_CPSO_H = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_CPSO_H(ii,1) = k2(end)+1;
    genra_CPSO_H(ii,2) = data_CPSO_H(ii,k2(end)+1);
end

X1 = genra_CPSO_H(:,1);
Y1 = genra_CPSO_H(:,2);
[Q1,~] = find(Y1 < Tao_1);
X1_ind = X1(Q1);
val_1 = Y1(Q1);



%% 2
diff_CPSO_S = diff(data_CPSO_S',1);
diff_CPSO_S(find(abs(diff_CPSO_S) < limit_de)) = 0;
diff_CPSO_S = diff_CPSO_S';
genra_CPSO_S = [];
for ii = 1 : size(diff_CPSO_S,1)
    temp_now = diff_CPSO_S(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_CPSO_S(ii,1) = k2(end)+1;
    genra_CPSO_S(ii,2) = data_CPSO_S(ii,k2(end)+1);
end
X2 = genra_CPSO_S(:,1);
Y2 = genra_CPSO_S(:,2);
[Q2,~] = find(Y2 < Tao_1);
X2_ind = X2(Q2);
val_2 = Y2(Q2);

   
%% 3
diff_DPSO = diff(data_DPSO',1);
diff_DPSO(find(abs(diff_DPSO) < limit_de)) = 0;
diff_DPSO = diff_DPSO';
genra_DPSO = [];
for ii = 1 : size(diff_DPSO,1)
    temp_now = diff_DPSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_DPSO(ii,1) = k2(end)+1;
    genra_DPSO(ii,2) = data_DPSO(ii,k2(end)+1);
end
X3 = genra_DPSO(:,1);
Y3 = genra_DPSO(:,2);
[Q3,~] = find(Y3 < Tao_1);
X3_ind = X3(Q3);
val_3 = Y3(Q3);


%% 4
diff_MPSO = diff(data_MPSO',1);
diff_MPSO(find(abs(diff_MPSO) < limit_de)) = 0;
diff_MPSO = diff_MPSO';
genra_MPSO = [];
for ii = 1 : size(diff_MPSO,1)
    temp_now = diff_MPSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_MPSO(ii,1) = k2(end)+1;
    genra_MPSO(ii,2) = data_MPSO(ii,k2(end)+1);
end
X4 = genra_MPSO(:,1);
Y4 = genra_MPSO(:,2);
[Q4,~] = find(Y4 < Tao_1);
X4_ind = X4(Q4);
val_4 = Y4(Q4);
if length(Q4)==0
    X4_ind =2;
    val_4 = min(Y4);
end


%% 5
diff_PSO = diff(data_PSO',1);
diff_PSO(find(abs(diff_PSO) < limit_de)) = 0;
diff_PSO = diff_PSO';
genra_PSO = [];
for ii = 1 : size(diff_PSO,1)
    temp_now = diff_PSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_PSO(ii,1) = k2(end)+1;
    genra_PSO(ii,2) = data_PSO(ii,k2(end)+1);
end
X5 = genra_PSO(:,1);
Y5 = genra_PSO(:,2);
[Q5,~] = find(Y5 < Tao_1);
X5_ind = X5(Q5);
val_5 = Y5(Q5);
if length(Q5)==0
    X5_ind =2;
    val_5 = min(Y5);
end
%% 6
diff_CLPSO = diff(data_CLPSO',1);
diff_CLPSO(find(abs(diff_CLPSO) < limit_de)) = 0;
diff_CLPSO = diff_CLPSO';
genra_CLPSO = [];
for ii = 1 : size(diff_CLPSO,1)
    temp_now = diff_CLPSO(ii,:);
    k6 = find(temp_now);
        if isempty(k6)
        k6(1) = 0;
    end
    genra_CLPSO(ii,1) = k6(end)+1;
    genra_CLPSO(ii,2) = data_CLPSO(ii,k6(end)+1);
end
X6 = genra_CLPSO(:,1);
Y6 = genra_CLPSO(:,2);
[Q6,~] = find(Y6 < Tao_1);
X6_ind = X6(Q6);
val_6 = Y6(Q6);
if length(Q6)==0
    X6_ind =2;
    val_6 = min(Y6)
end


%% 7
diff_proposed = diff(data_proposed',1);
diff_proposed(find(abs(diff_proposed) < limit_de)) = 0;
diff_proposed = diff_proposed';
genra_proposed = [];
for ii = 1 : size(diff_proposed,1)
    temp_now = diff_proposed(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_proposed(ii,1) = k2(end)+1;
    genra_proposed(ii,2) = data_proposed(ii,k2(end)+1);
end

X7 = genra_proposed(:,1);
Y7 = genra_proposed(:,2);
[Q7,~] = find(Y7 < Tao_1);
X7_ind = X7(Q7);
val_7 = Y7(Q7);


%% 8
diff_NPSO = diff(data_NPSO',1);
diff_NPSO(find(abs(diff_NPSO) < limit_de)) = 0;
diff_NPSO = diff_NPSO';
genra_NPSO = [];
for ii = 1 : size(diff_NPSO,1)
    temp_now = diff_NPSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_NPSO(ii,1) = k2(end)+1;
    genra_NPSO(ii,2) = data_NPSO(ii,k2(end)+1);
end

X8 = genra_NPSO(:,1);
Y8 = genra_NPSO(:,2);
[Q8,~] = find(Y8 < Tao_1);
X8_ind = X8(Q8);
val_8 = Y8(Q8);

%% 9
diff_GC_PSO = diff(data_GC_PSOd',1);
diff_GC_PSO(find(abs(diff_GC_PSO) < limit_de)) = 0;
diff_GC_PSO = diff_GC_PSO';
genra_NPSO = [];
for ii = 1 : size(diff_GC_PSO,1)
    temp_now = diff_GC_PSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_GC_PSOd(ii,1) = k2(end)+1;
    genra_GC_PSOd(ii,2) = data_GC_PSOd(ii,k2(end)+1);
end

X9 = genra_GC_PSOd(:,1);
Y9 = genra_GC_PSOd(:,2);
[Q9,~] = find(Y9 < Tao_1);
X9_ind = X9(Q9);
val_9 = Y9(Q9);
if isempty(Q9)
    X9_ind =1;
    val_9 = min(Y9);
end

%% 10
diff_DNSPSO = diff(data_DNSPSO',1);
diff_DNSPSO(find(abs(diff_DNSPSO) < limit_de)) = 0;
diff_DNSPSO = diff_DNSPSO';
genra_DNSPSO = [];
for ii = 1 : size(diff_DNSPSO,1)
    temp_now = diff_DNSPSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_DNSPSO(ii,1) = k2(end)+1;
    genra_DNSPSO(ii,2) = data_DNSPSO(ii,k2(end)+1);
end

X10 = genra_DNSPSO(:,1);
Y10 = genra_DNSPSO(:,2);
[Q10,~] = find(Y10 < Tao_1);
X10_ind = X10(Q10);
val_10 = Y10(Q10);
if isempty(Q10)
    X10_ind =1;
    val_10 = min(Y10);
end

%% 11
diff_MPCPSO = diff(data_MPCPSO',1);
diff_MPCPSO(find(abs(diff_MPCPSO) < limit_de)) = 0;
diff_MPCPSO = diff_MPCPSO';
genra_MPCPSO = [];
for ii = 1 : size(diff_MPCPSO,1)
    temp_now = diff_MPCPSO(ii,:);
    k2 = find(temp_now);
        if isempty(k2)
        k2(1) = 0;
    end
    genra_MPCPSO(ii,1) = k2(end)+1;
    genra_MPCPSO(ii,2) = data_MPCPSO(ii,k2(end)+1);
end
X11 = genra_MPCPSO(:,1);
Y11 = genra_MPCPSO(:,2);
[Q11,~] = find(Y11 < Tao_1);
X11_ind = X11(Q11);
val_11= Y11(Q11);
if isempty(Q11)
    X11_ind =1;
    val_11 = min(Y11);
end
%% 12
diff_LPSO = diff(data_LPSO',1);
diff_LPSO(find(abs(diff_LPSO) < limit_de)) = 0;
diff_LPSO = diff_LPSO';
genra_LPSO = [];
for ii = 1 : size(diff_LPSO,1)
    temp_now = diff_LPSO(ii,:);
    k12 = find(temp_now);
        if isempty(k12)
        k12(1) = 0;
    end
    genra_LPSO(ii,1) = k12(end)+1;
    genra_LPSO(ii,2) = data_LPSO(ii,k12(end)+1);
end
X12 = genra_LPSO(:,1);
Y12 = genra_LPSO(:,2);
[Q12,~] = find(Y12 < Tao_1);
X12_ind = X12(Q12);
val_12= Y12(Q12);
if isempty(Q12)
    X12_ind =1;
    val_12 = min(Y12);
end


%% 14
diff_LPSO_EVOLER = diff(data_LPSO_EVOLER',1);
diff_LPSO_EVOLER(find(abs(diff_LPSO_EVOLER) < limit_de)) = 0;
diff_LPSO_EVOLER = diff_LPSO_EVOLER';
genra_LPSO_EVOLER = [];
for ii = 1 : size(diff_LPSO_EVOLER,1)
    temp_now = diff_LPSO_EVOLER(ii,:);
    k14 = find(temp_now);
        if isempty(k14)
        k14(1) = 0;
    end
    genra_LPSO_EVOLER(ii,1) = k14(end)+1;
    genra_LPSO_EVOLER(ii,2) = data_LPSO_EVOLER(ii,k14(end)+1);
end
X14 = genra_LPSO_EVOLER(:,1);
Y14 = genra_LPSO_EVOLER(:,2);
[Q14,~] = find(Y14 < Tao_1);
X14_ind = X14(Q14);
val_14= Y14(Q14);
if isempty(Q14)
    X14_ind =1;
    val_14 = min(Y14);
end

%% 15
diff_downhill_EVOLER = diff(data_downhill_EVOLER',1);
diff_downhill_EVOLER(find(abs(diff_downhill_EVOLER) < limit_de)) = 0;
diff_downhill_EVOLER = diff_downhill_EVOLER';
genra_downhill_EVOLER = [];
for ii = 1 : size(diff_downhill_EVOLER,1)
    temp_now = diff_downhill_EVOLER(ii,:);
    k15 = find(temp_now);
        if isempty(k15)
        k15(1) = 0;
    end
    genra_downhill_EVOLER(ii,1) = k15(end)+1;
    genra_downhill_EVOLER(ii,2) = data_downhill_EVOLER(ii,k15(end)+1);
end
X15 = genra_downhill_EVOLER(:,1);
Y15 = genra_downhill_EVOLER(:,2);
[Q15,~] = find(Y15 < Tao_1);
X15_ind = X15(Q15);
val_15= Y15(Q15);
if isempty(Q15)
    X15_ind =1;
    val_15 = min(Y15);
end


Num = 50;
X_PSO = X5_ind;
[y_pso,x_pso]=hist(X_PSO,20);

X_DPSO = X3_ind;
[y_dpso,x_dpso]=hist(X_DPSO,20);

X_CPSOS = X2_ind;
[y_cpsos,x_cpsos]=hist(X_CPSOS,20);

X_CPSOH = X1_ind;
[y_cpsoh,x_cpsoh]=hist(X_CPSOH,20);

X_CLPSO = X6_ind;
[y_clpso,x_clpso]=hist(X_CLPSO,20);

X_MLPSO = X4_ind;
[y_mlpso,x_mlpso]=hist(X_MLPSO,20);

X_New_PSO1 = X7_ind;
[y_new_pso1,x_new_pso1]=hist(X_New_PSO1,20);


X_NPPSO = X8_ind;
[y_nppso,x_nppso]=hist(X_NPPSO,10);

X_GCPSO = X9_ind;
[y_gcpso,x_gcpso]=hist(X_GCPSO,10);

X_DNSPSO = X10_ind;
[y_dnspso,x_dnspso]=hist(X_DNSPSO,20);

X_MPCPSO = X11_ind;
[y_mpcpso,x_mpcpso]=hist(X_MPCPSO,20);

X_LPSO = X12_ind;
[y_lpso,x_lpso]=hist(X_LPSO,20);

X_LPSO_EVOLER = X14_ind;
[y_lpso_evoler,x_lpso_evoler]=hist(X_LPSO_EVOLER,20);

X_downhill_EVOLER = X15_ind;
[y_downhill_evoler,x_downhill_evoler]=hist(X_downhill_EVOLER,20);



%% Figure 4-g: achievable minimal costs of different methods, and the probabilities of finding them
figure
mean_various = [min(val_1), min(val_2), min(val_3), min(val_4) ,min(val_5), min(val_6), min(val_7), min(val_8), min(val_9), min(val_10), min(val_11), min(val_12), min(val_14), min(val_15)];
Probablity = [length(X_CPSOH)/Num length(X_CPSOS)/Num length(X_DPSO)/Num length(X_MLPSO)/Num length(X_PSO)/Num  length(X_CLPSO)/Num length(X_New_PSO1)/Num length(X_NPPSO)/Num length(X_GCPSO)/Num length(X_DNSPSO)/Num length(X_MPCPSO)/Num length(X_LPSO)/Num  length(X_LPSO_EVOLER)/Num length(X_downhill_EVOLER)/Num];

h14= plot(mean_various(14),Probablity(14),'s','MarkerSize',8,'lineWidth',1.5,'color',[0.9 0.4 0.2]);
set(h14,'MarkerFaceColor',get(h14,'color'));%downhill_EVOLER
hold on;
h13= plot(mean_various(13),Probablity(13),'*','MarkerSize',12 ,'lineWidth',1.5,'color',[0.1 0.1 0.7]);
set(h13,'MarkerFaceColor',get(h13,'color'));%LPSO_EVOLER
hold on;
h7 = plot(mean_various(7),Probablity(7),'pentagram','MarkerSize',8,'lineWidth',2,'color',[1.00 0.00 0.00])
set(h7,'MarkerFaceColor',get(h7,'color'));%EVOLER
hold on
h12= plot(mean_various(12),Probablity(12),'d','MarkerSize',8,'lineWidth',1.5,'color',[0.2 0.6 0.4])
set(h12,'MarkerFaceColor',get(h12,'color'));%LPSO
hold on;
h10 = plot(mean_various(10),Probablity(10),'<','MarkerSize',10,'lineWidth',1.5,'color',[0.72,0.27,1.00])
set(h10,'MarkerFaceColor',get(h10,'color'));%DNSPSO
hold on
h11= plot(mean_various(11),Probablity(11),'>','MarkerSize',10,'lineWidth',1.5,'color',[0.07,0.62,1.00])
set(h11,'MarkerFaceColor',get(h11,'color'));%MPCPSO
hold on;
h9 = plot(mean_various(9),Probablity(9),'hexagram','MarkerSize',10,'lineWidth',2,'color',[0.00,0.00,1.00])
set(h9,'MarkerFaceColor',get(h9,'color'));%GC_PSO
hold on
h8 = plot(mean_various(8),Probablity(8),'x','MarkerSize',10,'lineWidth',2,'color',[0.4 0.7 0.00])
set(h8,'MarkerFaceColor',get(h8,'color'));%NPSO
hold on
h4 = plot(mean_various(4),Probablity(4),'*','MarkerSize',14,'lineWidth',1.5,'color',[1.00,0.84,0.00])
set(h4,'MarkerFaceColor',get(h4,'color'));%MCPSO
hold on
h6 = plot(mean_various(6),Probablity(6),'*','MarkerSize',10,'lineWidth',1.5,'color',[1.00 0.27 0.01])
set(h6,'MarkerFaceColor',get(h6,'color'));%CLPSO
hold on
h1 = plot(mean_various(1),Probablity(1),'Square','MarkerSize',8,'lineWidth',1.5,'color',[0.20,0.80,0.20])
set(h1,'MarkerFaceColor',get(h1,'color'));%CPSO-H
hold on
h2 = plot(mean_various(2),Probablity(2),'diamond','MarkerSize',10,'lineWidth',1.5,'color',[0.50,1.00,0.00])
set(h2,'MarkerFaceColor',get(h2,'color'));%CPSO-S
hold on
h3 = plot(mean_various(3),Probablity(3),'x','MarkerSize',12,'lineWidth',1.5,'color',[0.00,0.78,0.55])
set(h3,'MarkerFaceColor',get(h3,'color'));%DPSO
hold on
h5 = plot(mean_various(5),Probablity(5),'+','MarkerSize',12,'lineWidth',1.5,'color',[1.00,0.27,0.01])
set(h5,'MarkerFaceColor',get(h5,'color'));%PSO
hold on
h1 = legend('EVOLER_d','EVOLER_l','EVOLER','Local PSO','DNSPSO','MPCPSO','GCPSO','NPSO','MCPSO','CLPSO','CPSO-H','CPSO-S','DPSO','PSO')

% box off
xlabel('Best Cost')
ylabel('Probablity of Global Solution')
% xlim([1.5302e4, 1.5410e4])

