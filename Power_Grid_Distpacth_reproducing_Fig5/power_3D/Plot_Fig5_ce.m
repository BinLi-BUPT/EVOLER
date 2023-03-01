% Copyright 2018, All Rights Reserved
% Code by Bin Li, Ziping Wei

clear all
clc
close all;

%% Running several files to generate the following data 
data_1 = load('CPSO_H_3D');              % run the MATLAB file 'CPSO_H_3D.m' to generate this data
data_2 = load('CPSO_S_3D');              % run the MATLAB file 'CPSO_S_3D.m' to generate this data
data_3 = load('DPSO_3D');                % run the MATLAB file 'DPSO_3D.m' to generate this data
data_4 = load('MCPSO_3D');               % run the MATLAB file 'MPSO_3D.m' to generate this data
data_5 = load('PSO_3D');                 % run the MATLAB file 'PSO_3D.m' to generate this data
data_6 = load('CLPSO_3D');               % run the MATLAB file 'CLPSO_3D.m' to generate this data
data_7 = load('EVOLER_3D');              % run the MATLAB file 'EVOLER_3D.m' to generate this data
data_8 = load('NPSO_3D');                % run the MATLAB file 'NPSO_3D.m' to generate this data
data_9 = load('GC_PSO_3D');              % run the MATLAB file 'GC_PSO_3D.m' to generate this data
data_10 = load('DNSPSO_3D');             % run the MATLAB file 'DNS_PSO_3D.m' to generate this data
data_11 = load('MPCPSO_3D');             % run the MATLAB file 'MPC_PSO_3D.m' to generate this data
data_12 = load('LPSO_3D');               % run the MATLAB file 'LPSO_3D.m' to generate this data
data_14 = load('LPSO_EVOLER_3D');        % run the MATLAB file 'EVOLER_3D.m' to generate this data
data_15 = load('downhill_EVOLER_3D');    % run the MATLAB file 'EVOLER_3D.m' to generate this data

data_CPSO_H = data_1.COSTS;
data_CPSO_S = data_2.COSTS;
data_DPSO = data_3.COSTS;
data_MPSO = data_4.COSTS;
data_PSO = data_5.COSTS;
data_CLPSO = data_6.COSTS;
data_proposed = data_7.Y;
data_NP_PSO = data_8.Y;
data_GC_PSOd = data_9.Y;
data_DNSPSO = data_10.Y;
data_MPCPSO = data_11.Y;
data_LPSO=data_12.COSTS;
data_LPSO_EVOLER=data_14.Y;
data_downhill_EVOLER=data_15.Y;

%==========================================================================
% figure 4-c: Averaged results of various evolutionary methods
K = 20;
N_add=2;
maxgen=150;
D=2;
marker_idx=1:20:(maxgen);

figure
x=1:maxgen;
x(1:N_add)=x(1:N_add)*K;
x(N_add+1:maxgen)=x(N_add)+(x(N_add+1:end)-N_add)*D;
y=mean(data_downhill_EVOLER);
y=y(1:maxgen);
semilogx(x,y,'-','LineWidth',2,'Color',[0.90,0.63,0.75],'MarkerIndices',marker_idx);
hold on;

x=1:maxgen;
x=x*K;
y=mean(data_LPSO_EVOLER);
y=y(1:maxgen);
semilogx(x,y,'--x','LineWidth',2,'Color',[1.00,0.28,0.75],'MarkerIndices',marker_idx);
hold on;

x=1:maxgen;
x=x*K;
y=mean(data_proposed);
y=y(1:maxgen);
semilogx(x,y,':ro','LineWidth',1.5,'MarkerIndices',marker_idx+5);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_LPSO);
y=y(1:maxgen);
semilogx(x,y,'-','LineWidth',1.5,'Color',[0.94,0.50,0.94],'Marker','pentagram','MarkerIndices',marker_idx);
hold on;

x=1:maxgen;
x=x*2*K;
y=mean(data_DNSPSO);
y=y(1:maxgen);
semilogx(x,y,':','LineWidth',2,'Color',[0.72,0.27,1.00],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*2*K;
y=mean(data_MPCPSO);
y=y(1:maxgen);
semilogx(x,y,'-.','LineWidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_GC_PSOd);
y=y(1:maxgen);
semilogx(x,y,'--','LineWidth',2,'Color',[0.22,0.22,0.90],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_NP_PSO);
y=y(1:maxgen);
semilogx(x,y,'-*','LineWidth',1.5,'Color',[1.00,0.45,0.79],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_MPSO);
y=y(1:maxgen);
semilogx(x,y,'-.','LineWidth',1.5,'Color',[0.93,0.69,0.13],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_CLPSO);
y=y(1:maxgen);
semilogx(x,y,'-','LineWidth',1,'Color',[0,0,0],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*1.5*K;
y=mean(data_CPSO_H);
y=y(1:maxgen);
semilogx(x,y,'-x','LineWidth',1,'Color',[0.25,0.78,0.25],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_CPSO_S);
y=y(1:maxgen);
semilogx(x,y,'-o','LineWidth',1,'Color',[0.52,1.00,0.17],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_DPSO);
y=y(1:maxgen);
semilogx(x,y,'-*','LineWidth',1,'Color',[0.02,0.72,0.58],'MarkerIndices',marker_idx);
hold on

x=1:maxgen;
x=x*K;
y=mean(data_PSO);
y=y(1:maxgen);
semilogx(x,y,'-.','LineWidth',2,'Color',[1.00,0.45,0.16],'MarkerIndices',marker_idx);
hold on

xlim([K,3000]);
ylim([8230,8480])
xlabel('Number of function evaluations')
ylabel('fitness')
% ylim([8230,8600])
h1 = legend('EVOLER_d','EVOLER_l','EVOLER','Local PSO','DNSPSO','MPCPSO','GCPSO','NPSO','MCPSO','CLPSO','CPSO-H','CPSO-S','DPSO','PSO')
set(legend,'EdgeColor',[1.00,1.00,1.00]);
limit_de = 1e-3;

%% 1
diff_CPSO_H = diff(data_CPSO_H',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_CPSO_H = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_CPSO_H(ii,1) = k2(end)+1;
    genra_CPSO_H(ii,2) = data_CPSO_H(ii,k2(end)+1);
end

X1 = genra_CPSO_H(:,1);
Y1 = genra_CPSO_H(:,2);



%% 2
diff_CPSO_H = diff(data_CPSO_S',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_CPSO_S = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_CPSO_S(ii,1) = k2(end)+1;
    genra_CPSO_S(ii,2) = data_CPSO_S(ii,k2(end)+1);
end
X2 = genra_CPSO_S(:,1);
Y2 = genra_CPSO_S(:,2);


   
%% 3
diff_CPSO_H = diff(data_DPSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_DPSO = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_DPSO(ii,1) = k2(end)+1;
    genra_DPSO(ii,2) = data_DPSO(ii,k2(end)+1);
end
X3 = genra_DPSO(:,1);
Y3 = genra_DPSO(:,2);


%% 4
diff_CPSO_H = diff(data_MPSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_MPSO = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_MPSO(ii,1) = k2(end)+1;
    genra_MPSO(ii,2) = data_MPSO(ii,k2(end)+1);
end
X4 = genra_MPSO(:,1);
Y4 = genra_MPSO(:,2);


%% 5
diff_CPSO_H = diff(data_PSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_PSO = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_PSO(ii,1) = k2(end)+1;
    genra_PSO(ii,2) = data_PSO(ii,k2(end)+1);
end
X5 = genra_PSO(:,1);
Y5 = genra_PSO(:,2);

%% 6
diff_CLPSO = diff(data_CLPSO',1);
diff_CLPSO(find(abs(diff_CLPSO) < limit_de)) = 0;
diff_CLPSO = diff_CLPSO';
genra_CLPSO = [];
for ii = 1 : size(diff_CLPSO,1)
    temp_now_clpso = diff_CLPSO(ii,:);
    k2 = find(temp_now_clpso);
    if length(k2)==0
        genra_CLPSO(ii,1)=1;
        genra_CLPSO(ii,2)=data_CLPSO(ii,1);
    else
        genra_CLPSO(ii,1) = k2(end)+1;
        genra_CLPSO(ii,2) = data_CLPSO(ii,k2(end)+1);
    end
end
% genra_CLPSO(find(genra_CLPSO(:,2) > 1e4),:) = [];
X6 = genra_CLPSO(:,1);
Y6 = genra_CLPSO(:,2);


%% 7
diff_CPSO_H = diff(data_proposed',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_proposed = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
    genra_proposed(ii,1) = k2(end)+1;
    genra_proposed(ii,2) = data_proposed(ii,k2(end)+1);
end

X7 = genra_proposed(:,1);
Y7 = genra_proposed(:,2);


%% 8
diff_CPSO_H = diff(data_NP_PSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_NPPSO = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
            if isempty(k2)
        k2(1) = 0;
    end
    genra_NPPSO(ii,1) = k2(end)+1;
    genra_NPPSO(ii,2) = data_NP_PSO(ii,k2(end)+1);
end

X8 = genra_NPPSO(:,1);
Y8 = genra_NPPSO(:,2);
% data_NP_PSO = data_8.Y;
% data_GC_PSOd = data_9.Y;

%% 9
diff_CPSO_H = diff(data_GC_PSOd',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_GC_PS = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
            if isempty(k2)
        k2(1) = 0;
    end
    genra_GC_PS(ii,1) = k2(end)+1;
    genra_GC_PS(ii,2) = data_GC_PSOd(ii,k2(end)+1);
end

X9 = genra_GC_PS(:,1);
Y9 = genra_GC_PS(:,2);


%% 10
diff_CPSO_H = diff(data_DNSPSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_GC_PS = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
            if isempty(k2)
        k2(1) = 0;
    end
    genra_DNS_PS(ii,1) = k2(end)+1;
    genra_DNS_PS(ii,2) = data_DNSPSO(ii,k2(end)+1);
end

X10 = genra_DNS_PS(:,1);
Y10 = genra_DNS_PS(:,2);

%% 11
diff_CPSO_H = diff(data_MPCPSO',1);
diff_CPSO_H(find(abs(diff_CPSO_H) < limit_de)) = 0;
diff_CPSO_H = diff_CPSO_H';
genra_GC_PS = [];
for ii = 1 : size(diff_CPSO_H,1)
    temp_now = diff_CPSO_H(ii,:);
    k2 = find(temp_now);
            if isempty(k2)
        k2(1) = 0;
    end
    genra_MPC_PS(ii,1) = k2(end)+1;
    genra_MPC_PS(ii,2) = data_MPCPSO(ii,k2(end)+1);
end

X11 = genra_MPC_PS(:,1);
Y11 = genra_MPC_PS(:,2);

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
    genra_LPSO_EVOLER(ii,2) = data_LPSO(ii,k14(end)+1);
end

X13 = genra_LPSO_EVOLER(:,1);
Y13 = genra_LPSO_EVOLER(:,2);

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

X14 = genra_downhill_EVOLER(:,1);
Y14 = genra_downhill_EVOLER(:,2);


%% BOX PLOT

color_def=[0.52,1.00,0.17;     0.25,0.78,0.25;     0.02,0.72,0.58;
           0.93,0.69,0.13;     1.00,0.45,0.16;     0.15,0.15,0.15;
           1.00,0.00,0.00;     1.00,0.45,0.79;     0.00,0.00,1.00;
           0.72,0.27,1.00;     0.07,0.62,1.00;     0.10,0.90,0.33;
           1.00,0.30,0.50;     0.00,0.15,0.43];
x_def{1}='CPSO-H';  color_def(1,:)=[0.25,0.78,0.25];
x_def{2}='CPSO-S';  color_def(2,:)=[0.52,1.00,0.17];
x_def{3}='DPSO';    color_def(3,:)=[0.02,0.72,0.58];
x_def{4}='MCPSO';   color_def(4,:)=[0.93,0.69,0.13];
x_def{5}='PSO';     color_def(5,:)=[1.00,0.45,0.16];
x_def{6}='CLPSO';   color_def(6,:)=[0,0,0];
x_def{7}='EVOLER';  color_def(7,:)=[1,0,0];
x_def{8}='NPSO';    color_def(8,:)=[1.00,0.45,0.79];
x_def{9}='GC-PSO';  color_def(9,:)=[0.22,0.22,0.90];
x_def{10}='DNSPSO'; color_def(10,:)=[0.72,0.27,1.00];
x_def{11}='MPCPSO'; color_def(11,:)=[0.07,0.62,1.00];
x_def{12}='Local PSO';color_def(12,:)=[0.94,0.50,0.94];
x_def{13}='EVOLER_l';color_def(13,:)=[1.00,0.28,0.75];
x_def{14}='EVOLER_d';color_def(14,:)=[0.90,0.63,0.75];

y_def{1}=X1;    y_def{2}=X2;    y_def{3}=X3;
y_def{4}=X4;    y_def{5}=X5;    y_def{6}=X6;
y_def{7}=X7;    y_def{8}=X8;    y_def{9}=X9;
y_def{10}=X10;    y_def{11}=X11;y_def{12}=X12;
y_def{13}=X13;    y_def{14}=X14;
x={};y={};
[a1,a2] = sort([mean(X1),mean(X2),mean(X3),mean(X4),mean(X5),mean(X6),mean(X7),mean(X8),mean(X9),mean(X10),mean(X11),mean(X12),mean(X13),mean(X14)])
for i=1:size(x_def,2)
    x{i}=x_def{a2(i)};
    colors(i,:)=color_def(a2(i),:);
    y{i}=y_def{a2(i)};
end
rng('default') % for reproducibility
%% Set Spacing
binWidth = 20;  % histogram bin widths
hgapGrp = .15;   % horizontal gap between pairs of boxplot/histograms (normalized)
hgap = 0.06;      % horizontal gap between boxplot and hist (normalized)
%% Compute histogram counts & edges
hcounts = cell(size(x,2),2); 
for i = 1:size(x,2)
    [hcounts{i,1}, hcounts{i,2}] = histcounts(y{i},'BinWidth',binWidth); 
end
maxCount = max([hcounts{:,1}]);
%% Plot boxplotsGroup()
fig = figure();
ax = axes(fig); 
hold(ax,'on')
% Convert y (mxn matrix) to 1xn cell array of mx1 vectors, required by boxplotWidths
 
xInterval = 1; %x-interval is always 1 with boxplot groups
normwidth = (1-hgapGrp-hgap)/2;     
boxplotWidth = xInterval*normwidth; 

% Plot colored boxplots
bph = boxplotGroup(ax,y,'Widths',boxplotWidth,'OutlierSize',3,'PrimaryLabels',x,'Colors',colors);
xlabel('Method');
ylabel('generations');
set(findobj(bph.boxplotGroup,'-property','LineWidth'), 'LineWidth', 1) % increase line widths
%% Add vertical histograms (patches) with matching colors
xCoordinate = 1:size(x,2);  %x-positions is always 1:n with boxplot groups
histX0 = xCoordinate + boxplotWidth/2 + hgap;    % histogram base
maxHeight = xInterval*normwidth;       % max histogram height
patchHandles = gobjects(1,size(x,2)); 
for i = 1:size(x,2)
    % Normalize heights 
    height = hcounts{i,1}/maxCount*maxHeight;
    % Compute x and y coordinates 
    xm = [zeros(1,numel(height)); repelem(height,2,1); zeros(2,numel(height))] + histX0(i);
    yidx = [0 0 1 1 0]' + (1:numel(height));
    ym = hcounts{i,2}(yidx);
    % Plot patches
    patchHandles(i) = patch(xm(:),ym(:),colors(i,:),'EdgeColor',colors(i,:),'LineWidth',1,'FaceAlpha',.45);
end