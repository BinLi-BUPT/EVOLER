clc;
clear all;

marker_idx=1:20:200;
t_add=21000*4.26/20;%single
t_add=t_add/3600;

cd EVOLER_D
data1=load('MEAN');
x1=1:30;
x1=x1*50;
x1=t_add+x1/3600;
x1=[0,t_add,x1];
y1=data1.MEAN;
y1=[30,30,y1];
cd ..\
h1=plot(x1,y1,'-','Color',[0.90,0.63,0.75],'Linewidth',1);

hold on;

cd EVOLER_L
data2=load('MEAN');
x2=1:30;
x2=x2*223.2;
x2=t_add+x2/3600;
x2=[0,t_add,x2];
y2=data2.MEAN;
y2=[30,30,y2];
cd ..\
h2=plot(x2,y2,'-.','Linewidth',1,'Color',[1.00,0.28,0.75],'MarkerIndices',marker_idx);

hold on;

cd EVOLER
data3=load('MEAN');
x3=1:20;
x3=x3*218;
x3=t_add+x3/3600;
x3=[0,t_add,x3];
y3=data3.MEAN;
y3=[30,30,y3];
cd ..\
h3=plot(x3,y3,'-s','Linewidth',1,'Color','r','MarkerIndices',marker_idx);

hold on;

cd LPSO
data4=load('MEAN');
x4=1:200;
x4=x4*223.2;
x4=x4/3600;
y4=data4.MEAN;
cd ..\
H4=plot(x4,y4,'-p','Linewidth',1,'Color',[0.94,0.50,0.94],'MarkerIndices',marker_idx);
hold on;

cd DNSPSO
data5=load('MEAN');
x5=1:200;
x5=x5*441;
x5=x5/3600;
y5=data5.MEAN;
cd ..\
H5=plot(x5,y5,'--','Linewidth',2,'Color',[0.72,0.27,1.00],'MarkerIndices',marker_idx);
hold on;

cd MPCPSO
data7=load('MEAN');
x7=1:200;
x7=x7*171;
x7=x7/3600;
y7=data7.MEAN;
cd ..\
H7=plot(x7,y7,'o-','Linewidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',marker_idx);
hold on;

cd GCPSO
data8=load('MEAN');
x8=1:200;
x8=x8*97;
x8=x8/3600;
y8=data8.MEAN;
cd ..\
H8=plot(x8,y8,'v-','Linewidth',1,'Color',[0.21,0.21,0.85],'MarkerIndices',marker_idx);
hold on;

cd NPSO
data9=load('MEAN');
x9=1:200;
x9=x9*98;
x9=x9/3600;
y9=data9.MEAN;
cd ..\
h9=plot(x9,y9,'-.','Linewidth',1.5,'Color',[0.96,0.46,0.77],'MarkerIndices',marker_idx);
hold on;

cd MCPSO
data10=load('MEAN');
x10=1:200;
x10=x10*198;
x10=x10/3600;
y10=data10.MEAN;
cd ..\
h10=plot(x10,y10,'*-','Linewidth',1,'Color',[0.93,0.69,0.13],'MarkerIndices',marker_idx);
hold on;

cd CLPSO
data11=load('MEAN');
x11=1:200;
x11=x11*214;
x11=x11/3600;
y11=data11.MEAN;
cd ..\
h11=plot(x11,y11,'+-','Linewidth',1,'Color',[0.50,0.50,0.50],'MarkerIndices',marker_idx);
hold on;

cd CPSO-H
data12=load('MEAN');
x12=1:200;
x12=x12*349;
x12=x12/3600;
y12=data12.MEAN;
cd ..\
h12=plot(x12,y12,'x-','Linewidth',1,'Color',[0.57,0.92,0.57],'MarkerIndices',marker_idx);
hold on;

cd CPSO-S
data13=load('MEAN');
x13=1:200;
x13=x13*416;
x13=x13/3600;
y13=data13.MEAN;
cd ..\
h13=plot(x13,y13,'h-','Linewidth',1,'Color',[0.06,0.97,0.60],'MarkerIndices',marker_idx);
hold on;

cd DPSO
data14=load('MEAN');
x14=1:200;
x14=x14*228;
x14=x14/3600;
y14=data14.MEAN;
cd ..\
h14=plot(x14,y14,'v-','Linewidth',1,'Color',[0.08,0.71,0.68],'MarkerIndices',marker_idx);
hold on;

cd PSO
data15=load('MEAN');
x15=1:200;
x15=x15*218;
x15=x15/3600;
y15=data15.MEAN;
cd ..\
h15=plot(x15,y15,'d-','Linewidth',1,'Color',[1.00,0.41,0.16],'MarkerIndices',marker_idx);
hold on;

xlabel('Time Complexity(hours)');
ylabel('Fitness');
legend('EVOLER_d','EVOLER_l','EVOLER','Local PSO', 'DNSPSO', 'MPCPSO',' GCPSO', 'NPSO', 'MCPSO', 'CLPSO', 'CPSO-H', 'CPSO-S', 'DPSO', 'PSO');