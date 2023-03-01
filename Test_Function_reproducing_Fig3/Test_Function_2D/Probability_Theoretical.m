se1 = [18	20	9.00E-25	16	8
50	200	7.00E-22	4	3
11	20	5.00E-27	3	2
2	80	3.00E-23	3	2
3	70	9.00E-24	2	2
25	16	6.00E-17	3	3
18	20	4.00E-07	19	9
50	100	6.00E-23	7	6
25	3	1.00E-29	2	2
1.3	120	3.00E-23	7	6
1	50	2.00E-06	40	14
160	20	2.00E-24	2	2
220	20	2.00E-22	2	2
1	20	5.00E-29	2	2
1.3	3	9.00E-28	2	2
21	11	1.00E-27	2	2
22	3	2.00E-24	2	2
1.4	3	1.00E-30	2	2
2.2 3   1e+03       3   2    
4.6	20	1.00E-26	4	3   
];

delta_f = se1(:,1); 
M=se1(:,2);
de=se1(:,3); 
s_all = se1(:,4); 
r = se1(:,5); 
epsilon = 0.02 * r./exp(s_all./r); 
de_s=(M-s_all).*de;
f=(1-epsilon).*2/pi.*atan(delta_f./de_s); 
figure
plot(repmat(1,1,21),'r-o')
hold on
plot(f,'b-*')
ylim([0.5,1])
set(gca,'xTick',[1:21]);
set(gca,'xTickLabel',{'Ackley','Griewank','Levy','Non-continuous Rastrigin','Rastrigin','Rosenbrock','Rotated Ackley','Rotated Griewank','Rotated Non-continuous Rastrigin','Rotated Rastrigin','Rotated Weierstrass','Schwefel','Sphere','Weierstrass','Shifted Levy','Shifted Rastrigin','Shifted Sphere','Shifted weierstrass','Hybrid Composition 1','Hybrid Composition 2' });


