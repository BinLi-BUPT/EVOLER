clc
clear all;

cd EVOLER_L_1
data =  load('pso-result30');
results(1,:)=data.result(1:30);
fitnessgbest(1)=data.result(30);
cd ..\
cd EVOLER_L_2
data =  load('pso-result30');
results(2,:)=data.result(1:30);
fitnessgbest(2)=data.result(30);
cd ..\
cd EVOLER_L_3
data =  load('pso-result30');
results(3,:)=data.result(1:30);
fitnessgbest(3)=data.result(30);
cd ..\
cd EVOLER_L_4
data =  load('pso-result30');
results(4,:)=data.result(1:30);
fitnessgbest(4)=data.result(30);
cd ..\
cd EVOLER_L_5
data =  load('pso-result30');
results(5,:)=data.result(1:30);
fitnessgbest(5)=data.result(30);
cd ..\
for i=1:size(results,2)
    AVE(i)=mean(results(:,i));
end
figure(1);
plot(AVE);
xlabel('generations');
ylabel('mean fitness');

[~,ind]=sort(fitnessgbest);
mid=ind(3)
mid_fitnessgbest=fitnessgbest(mid)
cd EVOLER_L_1
t = load('data_double_peak.txt');
sd = load('pso-result30.mat');
for ii=1:30
    as = load(['pso-result',num2str(ii)]);
    trans_result = as.trans_result;
    trans_result2 = trans_result;
    target=t;
    for jj=1:50
        fitness(jj)=sum(abs(target-trans_result2(jj,:)));
    end
    [m,p]=min(fitness);
    trans_result_best(ii,:)=trans_result2(p,:);
    trans_result2(p,:);
    s_1(ii) = m;
end
[inde_1, inde_2] = min(s_1);
XData=linspace(4.0000e-07,1.0000e-06,201);
XData=XData*1e9;
figure(2);
plot(XData,trans_result_best(inde_2,:),'-+','LineWidth',1.4);
hold on;
plot(XData,target,'--k');
hold on
h = fill(XData,target,'r');
set(h,'edgealpha',0,'facealpha',0.1)

legend('designed','target');
xlabel('wavelength(nm)');
ylabel('Transmission');
cd ..\