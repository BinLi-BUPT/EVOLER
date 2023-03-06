% Copyright 2022, All Rights Reserved
% Code by Bin Li & Ziping Wei

clc
clear  all
clc
close all
rf2 = @(x)schewfel(x);    %[-500,500]


d = 2;                    % Problem Dimension
F_min =[];
F_Pattern=[];
F_SA=[];
F_GA=[];
F_PSO=[];
F_SUR = [];
F_Global=[];
% The probability of finding the global optimum with different determinstic and heuristic methods£º GA: genetic
% algorithm, PS: pattern search, SA: simulated annealing, SQP: sequential quadratic programming, Multi-start: SQP with multiple starts.
Num = 500;                                    % Number of simulation trials
for iii=1:500
    iii
    x_bound = 500;
    x0 = 2*x_bound*[rand(1,d)-0.5];           % random start point
    n = length(x0);                           % number of variables
    lb = [-x_bound * ones(1,n)];              % lower-bound of variables
    ub = [x_bound * ones(1,n)];               % upper-bound of variables
    %----------------------------------------------------------------------
    % Deterministic optimization method: SQP
    options = optimoptions(@fmincon,'Algorithm','sqp');
    [xf,ff,flf,of] = fmincon(rf2,x0,[],[],[],[],lb,ub,[],options);
    F_min=[F_min,ff];
    %----------------------------------------------------------------------
    % Deterministic optimization method: Pattern Search
    [xp,fp,flp,op] = patternsearch(rf2,x0,[],[],[],[],lb,ub);
    F_Pattern=[F_Pattern,fp];
    %----------------------------------------------------------------------
    % Heuristic optimization method: Simulated Annealing
    [x,fval,exitFlag,output] = simulannealbnd(rf2,x0,lb,ub);
    F_SA=[F_SA,fval];
    %----------------------------------------------------------------------
    % Heuristic optimization method: Genetic Algorithm
    initpop = 2*x_bound*[rand(50,n)-0.5];
    opts = optimoptions('ga','InitialPopulationMatrix',initpop);
    [xga,fga,flga,oga] = ga(rf2,n,[],[],[],[],lb,ub,[]);
    F_GA=[F_GA,fga];
    %----------------------------------------------------------------------
    % Heuristic optimization method: Particle Swarm Optimization
    opts = optimoptions('particleswarm','SwarmSize',50);
    [xpso,fpso,flgpso,opso] = particleswarm(rf2,n,lb,ub,opts);
    F_PSO=[F_PSO,fpso];
    %----------------------------------------------------------------------
    % Deterministic global optimization method: Surrogate Optimization
    opts = optimoptions('surrogateopt','PlotFcn',[],'MaxFunctionEvaluations',400);
    %  default MaxFunctionEvaluations=max(200,50*nvar)
    [xsur,fsur,flgsur,osur] = surrogateopt(rf2,lb,ub,opts);
    F_SUR=[F_SUR,fsur];
    %----------------------------------------------------------------------
    % Deterministic global optimization method: Multi-starts + SQP
    problem = createOptimProblem('fmincon','objective',rf2,...
        'x0',x0,'lb',lb,'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    gs = GlobalSearch;
    [xg,fg,flg,og] = run(gs,problem);
    F_Global=[F_Global,fg];
end


threshold_success = 3e-5; %-1200;%  1e-20; % Threshold 
p_min = length(find(F_min<threshold_success))/length(F_min)
p_Pattern = length(find(F_Pattern<threshold_success))/length(F_Pattern)
p_SA = length(find(F_SA<threshold_success))/length(F_SA)
p_GA = length(find(F_GA<threshold_success))/length(F_GA)
p_PSO = length(find(F_PSO<threshold_success))/length(F_PSO)
p_SUR = length(find(F_SUR<threshold_success))/length(F_SUR)
p_Global = length(find(F_Global<threshold_success))/length(F_Global)
P_steer = 1;

%------Fig2-g: The probability of different methods to find global optimum--------%
figure
a=[P_steer;p_Global;p_PSO;p_GA;p_SUR;p_Pattern;p_min;p_SA];
b=diag(a);
c=bar(b,'stack','BarWidth',0.6);
color=[0.92,0.51,0.75;
    0.00,0.00,1.00;
    0.46,0.46,0.98;
    0.30,0.75,0.93;
    0.41,0.67,0.77;
    0.56,0.69,0.74;
    0.52,0.52,0.58;
    0.45,0.45,0.47];
alpha_face = [0.8 0.5 0.5 0.5 0.5 0.5 0.5 0.8];
for i=1:8
    set(c(i),'FaceColor',color(i,:));
    set(c(i),'FaceAlpha',alpha_face(i));
    set(c(i),'EdgeColor','none');
    %alpha(0.5)
end
hold on;
plot([-2:1:9],ones(1,12),'r--');
set(gca,'xticklabel',[]);
box off;
xlim([0.3 length(a)+1-.3]);
ylim([0 1.1]);
ylabel('Probability of Global Optimum')
save('Global_optimization_comparison_schw_K50_globalSQP_finnal.mat')
