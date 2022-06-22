function Bed_Optimization(m,iter)
%Bed_Optimization: Find iter optimal beds and display results as heatmaps
%                  distributions using three methods
%                  Crude, Simulated annealing and Pattern search
%                  
%
% Input:
%   m               Total bed capacity in the hospital 
%
%   iter            number of the optimal solutions to be found


%% MATLAB Optimization Toolbox Patternsearch

PSoptions = optimoptions('patternsearch','Display','iter');%display settings

Objfcn = @(x) f(round(x(1)),round(x(2)),m); %defining objective function

%%%%Finding solution%%%%%%%
for i =1:iter
cA0 = randi(m); %random initial guess of bed capacity for ward A
X0 = [cA0,randi(m-cA0)]; %random initial guess of bed capacities
%find minimum using solver
[Xps(i,:),Fps] = patternsearch(Objfcn,X0,[],[],[],[],[1,1],[m,m],PSoptions);
end

%%%%%%%PLOT solutions%%%%%%%%%%%%%
figure;

histogram2(Xps(:,1),Xps(:,2),'BinMethod','integers','DisplayStyle','tile','ShowEmptyBins','on')
title("Pattern search optimal solution - "+string(m)+" beds")
xlabel("c_A")
ylabel("c_C")
colorbar
saveas(gcf,"PS_"+string(m)+"_beds.png")
saveas(gcf,"PS_"+string(m)+"_beds")
%% Brute force method


S_10 = zeros(m,m,10); %Matrix for simulation results

for i = 1:iter
for capA = 1:m
    for capC =1:m
        if capA+capC <=m %if total capacity is not exceeded
            S_10(capA,capC,i)=f(capA,capC,m); %simulate for the bed distribution
        else
            S_10(capA,capC,i) =NaN; %The capacity is exceeded - no simulation
        end
    end
end
i
end

%caluculate mean of iter simulations
S = mean(S_10,3);
% save results
load Optimiziation_75_beds.mat
%find the capacites corresponding to minimum average of iter simulations
[Cap_A,Cap_C] = find(S==min(S,[],'all'));
Cap_B=m-Cap_A-Cap_C;% caluclate capacity B

%plot average of rellocated patients
figure;
imagesc(S)
ax = gca;
ax.YDir = 'normal';
xlabel("c_A")
ylabel("c_C")
colorbar
saveas(gcf,"opt_heat_map_75beds_lognorm.png")

%plot contour plot of average of relocated patients
figure;
hold on
contour(S,'ShowText','on')
plot(Cap_A,Cap_C,'*r')
for i = 1:10
[Cap_A,Cap_C] = find(S_10(:,:,i)==min(S_10(:,:,i),[],'all'));
plot(Cap_A,Cap_C,'*b')
end

xlabel("c_A")
ylabel("c_C")
saveas(gcf,"opt_contour_map_75beds_lognorm.png")
hold off

save Optimiziation_75_beds.mat


%% Simulated annealing

T = @(k) 1/sqrt(1+k); %Cooling scheme
X_opt = zeros(iter,2); %Allocate space for optimal solutions


for j = 1:iter %simulate iter times

cA0 = randi(m); %random initial guess of bed capacity for ward A

X = [cA0,randi(m-cA0)]; %random initial guess of bed capacities

i=0; %iteration counter
flag = true; %convergence criterion flag
while flag; %run until converges
    i =i+1; %increase iteration counter
    Y = [100,100]; %ensure that random walk is initiated
    %%%%%Random walk%%%%%%
    dirct = mod(i,2); %systematicly switch directions
    while sum(Y)>m || any(Y<0) %run until allowed direction is chosen
         delta_X = [sign(randn)*dirct,sign(randn)*(1-dirct)]; %choose direction
         Y = X(i,:)+delta_X; %New candidate
    end
    %%%%%end of random walk%%%
    UY = f(Y(1),Y(2),m); %value of energy/cost function at candidate
    UX = f(X(i,1),X(i,2),m); %value of energy/cost function at current position

    if UY<UX
        X(i+1,:) = Y; %accept 100%
    else
        if rand()<= exp(-(UY-UX)/T(i)) %accept with probability
            X(i+1,:) = Y; 
        else
            X(i+1,:) = X(i,:); %reject
        end
    end

%%%%Check for convergence
if i>20
    %terminate if the position was not changed for last 10 iterations
    flag = sum(abs(diff(X(end-10:end,:))),'all')>0;
end
end
%update current position
X_opt(j,:) = X(end,:);
j %display progress
end

%Plot optimal solutions
figure;
histogram2(X_opt(:,1),X_opt(:,2),'BinMethod','integers','DisplayStyle','tile','ShowEmptyBins','on')
xlabel("c_A")
ylabel("c_C")
colorbar
title("Simulated annealing optimal solution - " +string(m)+" beds")
saveas(gcf,"SA_"+string(m)+"_beds.png")
saveas(gcf,"SA_"+string(m)+"_beds")

save("Opt_"+string(m)+"beds.mat")
end


function S = f(capA,capC,m)
% Objective function

% ward A parameters
mu1 = log(4*sqrt(2));
s2_1 = log(2); 


% ward B parameters
mu2 = log(6*sqrt(2));
s2_2 = log(2); 


% ward C parameters
mu3 = log(5*sqrt(2));
s2_3 = log(2);

%Check if the total capacity is not exceded
% Needed for pattern search
if capA + capC >m
    S = inf;
else

%simulate
[Rejec, Realloc, ~, ~] = BedUtil([capA,m-capA-capC,capC],[mu1,mu2,mu3],[s2_1, s2_2, s2_3]);

%calculate the cost
S = sum(Rejec,'all')+sum(Realloc,'all');
end
end