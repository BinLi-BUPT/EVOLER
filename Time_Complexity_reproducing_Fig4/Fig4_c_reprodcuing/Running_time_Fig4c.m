clear all
load('time_end_all.mat')
Time_all = [mean(time1),mean(time2),mean(time3),mean(time4),mean(time5),mean(time6),mean(time7),mean(time8),mean(time9),...
    mean(time10),mean(time11),mean(time_evoler),mean(time13),mean(time14)];
[a1,a2] = sort(Time_all);
Time_end = Time_all(a2);
Time_all = Time_end;
for i=1:14
    dataT(i,i)=Time_all(i);
end
data1=dataT(1,:);
data2=dataT(2,:);
data3=dataT(3,:);
data4=dataT(4,:);
data5=dataT(5,:);
data6=dataT(6,:);
data7=dataT(7,:);
data8=dataT(8,:);
data9=dataT(9,:);
data10=dataT(10,:);
data11=dataT(11,:);
data12=dataT(12,:);
data13=dataT(13,:);
data14=dataT(14,:);

figure;
set (gcf,'WindowStyle','normal','Position', [500,400,700,300],'color','w');
bar(data1,'FaceColor',[0.90,0.63,0.75],'BarWidth',0.8);
hold on
bar(data2,'FaceColor',[1,0,0],'BarWidth',0.8);
hold on
bar(data3,'FaceColor',[1.00,0.28,0.75],'BarWidth',0.8);
hold on
bar(data4,'FaceColor',[0.96,0.46,0.77],'BarWidth',0.8);
hold on
bar(data5,'FaceColor',[0.50,0.50,0.50],'BarWidth',0.8);
hold on
bar(data6,'FaceColor',[0.93,0.69,0.13],'BarWidth',0.8);
hold on
bar(data7,'FaceColor',[1.00,0.41,0.16],'BarWidth',0.8); 
hold on
bar(data8,'FaceColor',[0,0,1],'BarWidth',0.8);
hold on
bar(data9,'FaceColor',[0.08,0.71,0.68],'BarWidth',0.8);
hold on
bar(data10,'FaceColor',[0.21,0.21,0.85],'BarWidth',0.8);
hold on
bar(data11,'FaceColor',[0.07,0.62,1.00],'BarWidth',0.8);
hold on
bar(data12,'FaceColor',[0.57,0.92,0.57],'BarWidth',0.8);
hold on
bar(data13,'FaceColor',[0.06,0.97,0.60],'BarWidth',0.8);
hold on
bar(data14,'FaceColor',[0.72,0.27,1.00],'BarWidth',0.8);
ylabel('Time/s')
set(gca,'xTick',[1:14]);
set(gca,'xTickLabel',{'EVOLER_d','EVOLER','EVOLER_l','NPSO','CLPSO','MCPSO', 'PSO','DPSO','Local PSO','GCPSO','MPCPSO','CPSO-S','CPSO-H','DNSPSO'});
xlabel('Method')

Time_11 = [time14;time_evoler;time13;time7;time5;time6;time1;time2;time11;time8;time9;time3;time4;time10];
for i=1:14
    ti=Time_11(i,:);
    MEAN(i)=mean(ti);
    SD(i)=std(ti);
    SEM(i)=SD(i)/sqrt(length(ti));
end
[max_a,~] = max(Time_11,[],2);
[min_b,~] = min(Time_11,[],2);
lu = max_a - Time_all';
up= Time_all' - min_b;
hold on
x = 1:14;
% er = errorbar(x,Time_all',errlow,errhigh);    
er = errorbar(x,MEAN,SEM,SEM); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
box off



% MEAN
% SEM


