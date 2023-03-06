clear all
close all

load('ff_all.mat');
load('time_all.mat');
% load('CLPSO_data.mat');
Time_sec = 1.35;
time_mean = [mean(time1),mean(time2),mean(time3),mean(time4),mean(time5),mean(time6),mean(time7),mean(time8),mean(time9),mean(time10),mean(time11)]/1200;
N_add = 632;

yc1 = ff_end1;
xc1 = time_mean(1).*[1:1:size(yc1,2)];
figure
yc1_end = yc1(:,round(Time_sec/time_mean(1)));
hist(yc1_end,180)
h = findobj(gca,'Type','patch');
h.FaceColor = [1.00,0.41,0.16];
h.EdgeColor = [1.00,0.41,0.16];

figure
yc2 = ff_end2;
xc2 = time_mean(2).*[1:1:size(yc2,2)];
yc2_end = yc2(:,round(Time_sec/time_mean(2)));
% histogram(yc2_end,80,'FaceColor',[0.08,0.71,0.68],'EdgeColor',[0.08,0.71,0.68])
hist(yc2_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.08,0.71,0.68];
h.EdgeColor = [0.08,0.71,0.68];
% shadedErrorBar(xc2,yc2,[0.08,0.71,0.68],[0.8,0.21,0.3],{@mean,@std});

figure
yc3 = ff_end3;
xc3 = time_mean(3).*[1:1:size(yc3,2)];
yc3_end = yc3(:,end);
% histogram(yc3_end,80,'FaceColor',[0.57,0.92,0.57],'EdgeColor',[0.57,0.92,0.57])
hist(yc3_end,20)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.57,0.92,0.57];
h.EdgeColor = [0.57,0.92,0.57];
% shadedErrorBar(xc3,yc3,[0.57,0.92,0.57],[0.8,0.21,0.3],{@mean,@std}); 

figure
yc4 = ff_end4;
xc4 = time_mean(4).*[1:1:size(yc4,2)];
yc4_end = yc4(:,round(Time_sec/time_mean(4)));
% histogram(yc4_end,80,'FaceColor',[0.06,0.97,0.60],'EdgeColor',[0.06,0.97,0.60])
hist(yc4_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.06,0.97,0.60];
h.EdgeColor = [0.06,0.97,0.60];
% hold on
% shadedErrorBar(xc4,yc4,[0.06,0.97,0.60],[0.8,0.21,0.3],{@mean,@std});

figure
yc5 = ff_end5;
xc5 = time_mean(5).*[1:1:size(yc5,2)];
yc5_end = yc5(:,round(Time_sec/time_mean(5)));
% histogram(yc5_end,80,'FaceColor',[0.50,0.50,0.50],'EdgeColor',[0.50,0.50,0.50])
% hold on
% shadedErrorBar(xc5,yc5,[0.50,0.50,0.50],[0.8,0.21,0.3],{@mean,@std});
hist(yc5_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.50,0.50,0.50];
h.EdgeColor = [0.50,0.50,0.50];

figure
yc6 = ff_end6;
xc6 = time_mean(6).*[1:1:size(yc6,2)];
yc6_end = yc6(:,round(Time_sec/time_mean(6)));
% histogram(yc6_end,80,'FaceColor',[0.93,0.69,0.13],'EdgeColor',[0.93,0.69,0.13])
hist(yc6_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.93,0.69,0.13];
h.EdgeColor = [0.93,0.69,0.13];
% hold on
% shadedErrorBar(xc6,yc6,[0.93,0.69,0.13],[0.8,0.21,0.3],{@mean,@std});

figure
yc7 = ff_end7;
xc7 = time_mean(7).*[1:1:size(yc7,2)];
yc7_end = yc7(:,round(Time_sec/time_mean(7)));
% histogram(yc7_end,80,'FaceColor',[0.96,0.46,0.77],'EdgeColor',[0.96,0.46,0.77])
% hold on
% shadedErrorBar(xc7,yc7,[0.96,0.46,0.77],[0.8,0.21,0.3],{@mean,@std});
hist(yc7_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.96,0.46,0.77];
h.EdgeColor = [0.96,0.46,0.77];


figure
yc8 = ff_end8;
xc8 = time_mean(8).*[1:1:size(yc8,2)];
yc8_end = yc8(:,round(Time_sec/time_mean(8)));
% histogram(yc8_end,80,'FaceColor',[0.21,0.21,0.85],'EdgeColor',[0.21,0.21,0.85])
hist(yc8_end,200)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.21,0.21,0.85];
h.EdgeColor = [0.21,0.21,0.85];
% hold on
% shadedErrorBar(xc8,yc8,[0.21,0.21,0.85],[0.8,0.21,0.3],{@mean,@std});

figure
yc9 = ff_end10;
xc9 = time_mean(9).*[1:1:size(yc9,2)];
yc9_end = yc9(:,round(Time_sec/time_mean(9)));
% histogram(yc9_end,80,'FaceColor',[0.07,0.62,1.00],'EdgeColor',[0.07,0.62,1.00])
hist(yc9_end,100)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.07,0.62,1.00];
h.EdgeColor = [0.07,0.62,1.00];
% hold on
% shadedErrorBar(xc9,yc9,[0.07,0.62,1.00],[0.8,0.21,0.3],{@mean,@std});

figure
yc10 = ff_end11;
xc10 = time_mean(10).*[1:1:size(yc10,2)];
yc10_end = yc10(:,round(Time_sec/time_mean(10)));
% histogram(yc10_end,80,'FaceColor',[0.72,0.27,1.00],'EdgeColor',[0.72,0.27,1.00])
hist(yc10_end,100)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.72,0.27,1.00];
h.EdgeColor = [0.72,0.27,1.00];
% hold on
% shadedErrorBar(xc10,yc10,[0.72,0.27,1.00],[0.8,0.21,0.3],{@mean,@std});

figure
yc11 = ff_end12;
xc11 = time_mean(11).*[1:1:size(yc11,2)];
yc11_end = yc11(:,round(Time_sec/time_mean(11)));
% histogram(yc11_end,80,'FaceColor',[0.94,0.50,0.94],'EdgeColor',[0.94,0.50,0.94])
hist(yc11_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.94,0.50,0.94];
h.EdgeColor = [0.94,0.50,0.94];
% hold on
% shadedErrorBar(xc11,yc11,[0.94,0.50,0.94],[0.8,0.21,0.3],{@mean,@std});
figure
load('shifted_Pro_data_levy.mat')
yc12 = Fit_pro;
% xc12 = time_mean(12).*[1:1:size(yc12,2)];
load('time_evoler.mat');
xc12 = zeros(1,900);
xc12(1:N_add) = T1.*[1:N_add];
xc12(633:900) = T1*N_add + time_mean(1).*[1:900-N_add];
% xc12 = [, [633:900]*0.0012N_add43000000];

yc12_end = yc12(:,891);
hist(yc12_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [1,0,0];
h.EdgeColor = [1,0,0];
% hold on
% shadedErrorBar(xc12,yc12,[1,0,0],[0.8,0.21,0.3],{@mean,@std});

load('Local_data_levy.mat')
figure
yc13 = Fit_pro;
% xc13 = time_mean(13).*[1:1:size(yc13,2)];
xc13 = zeros(1,900);
xc13(1:N_add) = T1.*[1:N_add];
xc13(633:900) = T1*N_add + time_mean(11).*[1:900-N_add];

yc13_end = yc13(:,900);
hist(yc13_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [1.00,0.28,0.75];
h.EdgeColor = [1.00,0.28,0.75];
% hold on
% shadedErrorBar(xc13,yc13,[1.00,0.28,0.75],[0.8,0.21,0.3],{@mean,@std});

load('Downhill_Pro_data_levy.mat')
figure
yc14 = Fit_pro;
% xc14 = time_mean(14).*[1:1:size(yc14,2)];
load('time_evoler_d1.mat')
xc14 = zeros(1,900);
xc14(1:N_add) = T1.*[1:N_add];
xc14(633:900) = T1*N_add + T2.*[1:900-N_add];

yc14_end = yc14(:,900);
% histogram(yc14_end,80,'FaceColor',[0.90,0.63,0.75],'EdgeColor',[0.90,0.63,0.75])
hist(yc14_end,80)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.90,0.63,0.75];
h.EdgeColor = [0.90,0.63,0.75];
% hold on
% shadedErrorBar(xc14,yc14,[0.90,0.63,0.75],[0.8,0.21,0.3],{@mean,@std});


figure
semilogy(xc1,mean(yc1,1),'d-','Linewidth',1,'Color',[1.00,0.41,0.16],'MarkerIndices',[1:50:length(xc1)]);
hold on
semilogy(xc2,mean(yc2,1),'v-','Linewidth',1,'Color',[0.08,0.71,0.68],'MarkerIndices',[1:50:length(xc2)]);
hold on
semilogy(xc3,mean(yc3,1),'h-','Linewidth',1,'Color',[0.57,0.92,0.57],'MarkerIndices',[1:50:length(xc3)]);
hold on
semilogy(xc4,mean(yc4,1),'x-','Linewidth',1,'Color', [0.06,0.97,0.60],'MarkerIndices',[1:50:length(xc4)]);
hold on
semilogy(xc5,mean(yc5,1),'+-','Linewidth',1,'Color',[0.50,0.50,0.50],'MarkerIndices',[1:50:length(xc5)]);
hold on
semilogy(xc6,mean(yc6,1),'*-','Linewidth',1,'Color',[0.93,0.69,0.13],'MarkerIndices',[1:50:length(xc6)]);
hold on
semilogy(xc7,mean(yc7,1),'.-','Linewidth',1.5,'Color',[0.96,0.46,0.77],'MarkerIndices',[1:50:length(xc7)]);
hold on
semilogy(xc8,mean(yc8,1),'v-','Linewidth',1,'Color',[0.21,0.21,0.85],'MarkerIndices',[1:50:length(xc8)]);
hold on
semilogy(xc9,mean(yc9,1),'o-','Linewidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',[1:50:length(xc9)]);
hold on
semilogy(xc10,mean(yc10,1),'--','Linewidth',2,'Color',[0.72,0.27,1.00],'MarkerIndices',[1:50:length(xc10)]);
hold on
semilogy(xc11,mean(yc11,1),'-p','Linewidth',1,'Color',[0.94,0.50,0.94],'MarkerIndices',[1:50:length(xc11)]);
hold on
semilogy(xc12,mean(yc12,1),'-s','Linewidth',1,'Color','r','MarkerIndices',[1:50:length(xc12)]);
hold on
semilogy(xc13,mean(yc13,1),'-.','Linewidth',1,'Color',[1.00,0.28,0.75],'MarkerIndices',[1:50:length(xc13)]);
hold on
semilogy(xc14,mean(yc14,1),'-.','Linewidth',1,'Color',[0.90,0.63,0.75],'MarkerIndices',[1:50:length(xc14)]);
xlabel('Time/sec')
ylabel('Fitness')
legend('PSO','DPSO','CPSO-S','CPSO-H','CLPSO','MCPSO','NPSO','GCPSO','MPCPSO','DNSPSO','Local PSO','EVOLER','EVOLER_l','EVOLER_d','location','SouthWest');
xlim([0,Time_sec])
ylim([1e-18,4000])

hist_PSO = yc1_end;
hist_DPSO = yc2_end;
hist_CPSOS = yc3_end;
hist_CPSOH = yc4_end;
hist_CLPSO = yc5_end;
hist_MCPSO = yc6_end;
hist_NPSO = yc7_end;
hist_GCPSO = yc8_end;
hist_MPCPSO = yc9_end;
hist_DNSPSO = yc10_end;
hist_Local_PSO = yc11_end;
hist_EVOLER = yc12_end;
hist_EVOLER_L = yc13_end;
hist_EVIOLER_D = yc14_end;
save hist_all hist_PSO hist_DPSO hist_CPSOS hist_CPSOH hist_CLPSO hist_MCPSO hist_NPSO hist_GCPSO hist_MPCPSO hist_DNSPSO hist_Local_PSO hist_EVOLER hist_EVOLER_L hist_EVIOLER_D

PSO_data = ff_end1(:,1:1200);
DPSO_data = ff_end2(:,1:1200);
CPSOS_data = ff_end3(:,1:1200);
CPSOH_data = ff_end4(:,1:1200);
CLPSO_data = ff_end5(:,1:2600);
MCPSO_data = ff_end6(:,1:1200);
NPSO_data = ff_end7(:,1:1200);
GCPSO_data = ff_end8(:,1:1200);
MPCPSO_data = ff_end10(:,1:1200);
DNSPSO_data = ff_end11(:,1:1200);
Local_PSO_data = ff_end12(:,1:1200);
de1 = load('shifted_Pro_data_levy.mat');
de2 = load('Local_data_levy.mat');
de3 = load('Downhill_Pro_data_levy.mat');
EVOLER_data = de1.Fit_pro;
EVOLER_L_data = de2.Fit_pro;
EVOLER_D_data = de3.Fit_pro;
save data_all PSO_data DPSO_data CPSOS_data CPSOH_data CLPSO_data MCPSO_data NPSO_data GCPSO_data MPCPSO_data DNSPSO_data Local_PSO_data EVOLER_data EVOLER_L_data EVOLER_D_data
% 
