% clc;
clear all;
% close all;
c1 = 1.49;
c2 = 1.49;  
maxg = 1000;          
sizepop = 20;       
D= 3;                
popmax = 550;        
popmin = 50;    
Vmax = 20;
Vmin = -Vmax;
Total = 850;
P_i_min_all = [100,100,50];
P_i_max_all = [500,500,300];
P_1_text = P_i_min_all(1:D-1);
P_2_text = P_i_max_all(1:D-1);
DMin=P_i_min_all(end);
DMax=P_i_max_all(end);
pop = [];
pop2 = [];

for iii = 1 : 50
    iii
    for i=1:sizepop
        while(1)
            pop2(i,:)= popmax*rand(1,D-1);
            ind_1 = pop2(i,:) < P_1_text;
            ind_peak = find(ind_1 == 1);
            pop2(i,ind_peak) = P_1_text(ind_peak);
            
            ind_2 = pop2(i,:) > P_2_text;
            ind_peak2 = find(ind_2 == 1);
            pop2(i,ind_peak2) = P_2_text(ind_peak2);
            pop2(find(pop2 < 0)) = 1;
            if Total - sum(pop2(i,:))<=DMax && Total - sum(pop2(i,:))>=DMin, break;end
        end
        pop(i,:)= ([pop2(i,:), Total - sum(pop2(i,:))]);
        speed(i,:)=Vmax*rands(1,D-1);
        fitness(i)= quad_cf(pop(i,:));
    end
    pBest= pop2;
    fitnesspbest=fitness;
    [fitnessgbest bestindex]=min(fitness);
    gBest=(pop2(bestindex,:));



    %% GCPSO 8
    
    c1=1.49; c2=1.49; w=0.72;
    ft = 1;                                 %range control
    success = 15;                           %numbers of success
    failure = 5;                            %numbers of failure
    suNumber = 0;                           %count the number of success
    faNumber = 0;                           %count the number of failure
    speedKeep = zeros(1,D-1);               %keep the speed of best particle
    nMin = bestindex;
    xBest = pop2;                           %each particle's best position
    yBest = gBest;                          %best position for each particles in history
    xFit =  fitnesspbest;                   %each particle's best fittness
    yFit =  fitnessgbest;                   %best fittness in history
    for iIteration = 1:maxg
        for iParticles = 1:sizepop
            pop_exi = Total - sum(pop2(iParticles,:));
            pop_sum = [pop2(iParticles,:), pop_exi];
            if (pop_exi >  P_i_max_all(end)) || (pop_exi < P_i_min_all(end))
                xFitness = 1e10;
            else
                xFitness = quad_cf(pop_sum);
            end
            %current particle fittness
            if xFitness < xFit(iParticles)            % renew particle's best position
                xFit(iParticles) = xFitness;
                xBest(iParticles,:) = pop2(iParticles,:);
            end
        end
        if min(xFit)  < yFit                          % renew best fittness
            [yFit,nMin] = min(xFit);
            yBest = xBest(nMin,:);
            %%first insert
            suNumber = suNumber + 1;
            faNumber = 0;
        else
            faNumber = faNumber + 1;
            suNumber = 0;
        end
        
        %%iteration end,renew speed and position
        speed = (speed * w + c1 * rand * (xBest - pop2) + c2 * rand * (repmat(yBest,sizepop,1) - pop2));%renew speed
        speed(speed < Vmin ) = Vmin;
        speed(speed > Vmax) = Vmax;
        pop2 = pop2 + speed;          %renew position
        %third insert
        if suNumber > success
            ft = 2 * ft;
            failure = failure + 1;  %renew failure 
        elseif faNumber > failure
            ft = 0.5 * ft;
            success = success + 1;  %renew success 
        end
        pop2(nMin,:) = yBest(1,:) + w * speedKeep + ft .* (1 - 2 * rand(1,D-1));

        for j=1:sizepop
            ind_1 = pop2(j,:) < P_1_text;
            ind_peak = find(ind_1 == 1);
            pop2(j,ind_peak) = P_1_text(ind_peak);
            
            ind_2 = pop2(j,:) > P_2_text;
            ind_peak2 = find(ind_2 == 1);
            pop2(j,ind_peak2) = P_2_text(ind_peak2);
            pop2(find(pop2 < 0)) = 1;
        end
        result(iIteration) = yFit;
    end
    g_Best_end = [yBest, (Total - sum(yBest))];
    Y(iii,:) = result;
end
save('GC_PSO_3D','Y')

