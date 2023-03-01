clc
clear all;

cd EVOLER_D_1
data =  load('EVOLER_d-result30');
results(1,:)=data.result(1:30);
gbestfitness(1)=data.result(30);
cd ..\
cd EVOLER_D_2
data =  load('EVOLER_d-result30');
results(2,:)=data.result(1:30);
gbestfitness(2)=data.result(30);
cd ..\
cd EVOLER_D_3
data =  load('EVOLER_d-result30');
results(3,:)=data.result(1:30);
gbestfitness(3)=data.result(30);
cd ..\
cd EVOLER_D_4
data =  load('EVOLER_d-result30');
results(4,:)=data.result(1:30);
gbestfitness(4)=data.result(30);
cd ..\
cd EVOLER_D_5
data =  load('EVOLER_d-result30');
results(5,:)=data.result(1:30);
gbestfitness(5)=data.result(30);
cd ..\
for i=1:size(results,2)
    AVE(i)=mean(results(:,i));
end
figure;
plot(AVE);
xlabel('generations');
ylabel('mean fitness');
legend('downhill EVOLER')


[~,ind]=sort(gbestfitness);
mid=ind(3)
mid_fitnessgbest=gbestfitness(mid)
cd EVOLER_D_3
t = load('data_double_peak.txt');
sd = load('EVOLER_d-result30.mat');
for ii=1:30
    as = load(['EVOLER_d-result',num2str(ii)]);
    trans_result = as.trans_result;
    trans_result2 = trans_result;
    target=t;
    for jj=1:6
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