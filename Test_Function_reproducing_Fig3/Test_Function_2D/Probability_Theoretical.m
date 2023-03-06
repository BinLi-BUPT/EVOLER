close all
se1 = [15	20	4.91165049884953e-14	16	2;
50 200	1.36403148198268e-13	4	2;
2  20	7.16591917453260e-15	4	3;
60 3	0.000223856820195165	2	2;   
3  20	3.64795165236133e-15	3	2;
0.2 80	7.00947854226503e-14	3	2;
0.1 70	4.37189189473434e-14	2	2;  
5  20	5.16768246692093e-10	3	1;
15 20	0.000310413446927976	19	4;
0.2 200	1.10972355754648e-13	7	2;
25 3	1.22080433768218e-15	2	2;
1  120	4.82694782496912e-14	7	5;
1 50	5.13656453652931e-05	40	12;
120 20	8.41831770810180e-14	2	1;
35  3	5.02120842330479e-15	2	2;
1  11	4.66166928437654e-15	2	1;
1e3 3	9.23746879619637e-13	2	2;
2   3	1.40276526187604e-16	2	1;
220 20	7.75683138748747e-13	2	2;
15 20	3.63917244158412e-16	2	1;
];

delta_f = se1(:,1); 
M=se1(:,2);
de=se1(:,3); 
s_all = se1(:,4); 
r = se1(:,5); 
epsilon = 1 * r./exp(3.5*s_all./r);
% de_s=(M-s_all).*de;
de_s=de;
f=(1-epsilon).*2/pi.*atan(delta_f./de_s); 
% figure
% plot(repmat(1,1,20),'r-o')
% hold on
% plot(f,'b-*')
% ylim([0.5,1])
% set(gca,'xTick',[1:20]);
% set(gca,'xTickLabel',{'Ackley','Griewank','Hybrid Composition 1','Hybrid Composition 3' ,'Levy', 'Non-continuous Rastrigin','Rastrigin','Rosenbrock','Rotated Ackley','Rotated Griewank','Rotated Non-continuous Rastrigin','Rotated Rastrigin','Rotated Weierstrass','Schwefel','Shifted Levy','Shifted Rastrigin','Shifted Sphere','Shifted weierstrass','Spheref','Weierstrass'});

figure
prac=ones(1,20);
for i=1:20
    prob(i,:)=[f(i),prac(i)];
end
b=bar(prob);
ylim([0,1]);
ch = get(b,'children');
x=1:20;
set(gca,'XTick',x,'XTickLabel',{'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12','f13','f14','f15','f16','f17','f18','f19','f20'});
legend('Theoretical low-bound','Emperical results');
ylabel('Probability of convergence to global optimum');
box off;

% theo=[1,0.998176236,0.971789312,0.939603001,0.989504963,0.989504963,...
%     0.939605233,0.999972463,0.999986585,0.99999043,0.939605233,...
%     0.962767085,0.999864404,0.999088118,0.939605233,0.999088118,...
% 0.939605233,0.999088118,0.939605233,0.999088118]
% 
% 
% prac=ones(1,20);
% for i=1:20
%     prob(i,:)=[theo(i),prac(i)];
% end
% figure
% b=bar(prob);
% ylim([0.5,1]);
% ch = get(b,'children');
% set(gca,'XTickLabel',{'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10','f11','f12','f13','f14','f15','f16','f17','f18','f19','f20'});
% legend('Theoretical low-bound','Emperical results');
% ylabel('Probability of convergence to global optimum');
% box off;