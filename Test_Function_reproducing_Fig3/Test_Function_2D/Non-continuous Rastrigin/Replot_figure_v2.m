clc
h_line=get(gca,'Children');%get lfill(xdata,ydata,'p');inehandles
xdata=get(h_line,'Xdata');
ydata=get(h_line,'Ydata');
K = 50;

Ydata = zeros(14,500);
Xdata = zeros(14,500);
for iii=1:14
    Ydata(1:length( ydata{iii}),iii)= ydata{iii} ;
    Xdata(1:length( xdata{iii}),iii)= xdata{iii} ;
end 
marker_idx=1:30:length(Xdata(1:500,1))*K;




figure

%% EVOLER_Downhill
Xdata_de=Xdata(1:500,1);
N_add=6;
Xdata_de(1:N_add)=Xdata_de(1:N_add)*K;
Xdata_de(N_add+1:500)=Xdata_de(N_add)+(Xdata_de(N_add+1:end)-N_add)*2;
semilogy(Xdata_de,Ydata(1:500,1),'Color',[0.90,0.63,0.75],'Linewidth',1)
%legend('Ours')
hold on
%% EVOLER_local
semilogy(Xdata(1:500,2)*K,Ydata(1:500,2),'-.','Linewidth',1,'Color',[1.00,0.28,0.75],'MarkerIndices',marker_idx)
%legend('MPCPSO')

%% EVOLER_PSO
semilogy(Xdata(1:500,3)*K,Ydata(1:500,3),'-s','Linewidth',1,'Color','r','MarkerIndices',marker_idx)
%legend('NDPSO')

%% Local_PSO
semilogy(Xdata(1:500,4)*K,Ydata(1:500,4),'-p','Linewidth',1,'Color',[0.94,0.50,0.94],'MarkerIndices',marker_idx)
%legend('GCPSO')

%% DNSSO
semilogy(Xdata(1:500,5)*2*K,Ydata(1:500,5),'--','Linewidth',2,'Color',[0.72,0.27,1.00],'MarkerIndices',marker_idx)
%legend('NPSO')

%% MPCPSO
semilogy(Xdata(1:500,6)*2*K,Ydata(1:500,6),'o-','Linewidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',marker_idx)
%legend('MCPSO')

%% GCPSO
semilogy(Xdata(1:500,7)*K,Ydata(1:500,7),'v-','Linewidth',1,'Color',[0.21,0.21,0.85],'MarkerIndices',marker_idx)
%legend('CLPSO')

%% NPSO
semilogy(Xdata(1:500,8)*K,Ydata(1:500,8),'.-','Linewidth',1.5,'Color',[0.96,0.46,0.77],'MarkerIndices',marker_idx)
%legend('CPSO-H')

%% MCPSO
semilogy(Xdata(1:500,9)*K,Ydata(1:500,9),'*-','Linewidth',1,'Color',[0.93,0.69,0.13],'MarkerIndices',marker_idx)
%legend('CPSO-S')

%% CLPSO
semilogy(Xdata(1:500,10)*K,Ydata(1:500,10),'+-','Linewidth',1,'Color',[0.50,0.50,0.50],'MarkerIndices',marker_idx)
%legend('DPSO')

%% CPSO-H
semilogy(Xdata(1:500,11)*2*K,Ydata(1:500,11),'x-','Linewidth',1,'Color',[0.57,0.92,0.57],'MarkerIndices',marker_idx)

%% CPSO-S
semilogy(Xdata(1:500,12)*(2*K+K/2),Ydata(1:500,12),'h-','Linewidth',1,'Color',[0.06,0.97,0.60],'MarkerIndices',marker_idx)

%% DPSO
semilogy(Xdata(1:500,13)*K,Ydata(1:500,13),'v-','Linewidth',1,'Color',[0.08,0.71,0.68],'MarkerIndices',marker_idx)

%% PSO
semilogy(Xdata(1:500,14)*K,Ydata(1:500,14),'d-','Linewidth',1,'Color',[1.00,0.41,0.16],'MarkerIndices',marker_idx)
h1 = legend('EVOLER_d','EVOLER_l','EVOLER','Local PSO','DNSPSO','MPCPSO','GCPSO','NPSO','MCPSO','CLPSO','CPSO-H','CPSO-S','DPSO','PSO')
set(legend,'EdgeColor',[1.00,1.00,1.00]);
xlim([0,5e4]);
xlabel('Number of function evaluations');
ylabel('Fitness');
box off
ax = gca;
ax.XAxis.Exponent = 4;


%% =====================================
% h_line=get(gca,'Children');%get lfill(xdata,ydata,'p');inehandles
% xdata=get(h_line,'xdata');
% ydata=get(h_line,'ydata');
% marker_idx=1:50:length(xdata);
% plot(xdata,ydata,'-.','Linewidth',1,'Color',[0.07,0.62,1.00],'MarkerIndices',marker_idx)
