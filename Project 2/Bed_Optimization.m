%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utilization of hospital beds during epidemics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
close all


%% MATLAB Optimization Toolbox - doesn't work
% x = optimvar('x',Type='integer',LowerBound=1,UpperBound=75);
% y = optimvar('y', Type='integer',LowerBound=1,UpperBound=75);
% obj = f(x,y);
% prob = optimproblem('Objective',obj);
% cnstr = x+y<=75;
% prob.Constraints.constr = cnstr;
% x0.x = 25;
% x0.y = 25;
% show(prob)
% [sol,fval] = solve(prob,x0);
S = zeros(75,75,10);
%% Brute force method
for i = 1:10
for capA = 1:75
    for capC =1:75
        if capA+capC <=75
            S(capA,capC,i)=f(capA,capC);
        else
            S(capA,capC,i) =NaN;
        end
    end
end
i
end

%%
S_10 = S;
S = mean(S,3);
%%
load Optimiziation_75_beds.mat
[Cap_A,Cap_C] = find(S==min(S,[],'all'));
Cap_B=75-Cap_A-Cap_C;
figure;
imagesc(S)
ax = gca;
ax.YDir = 'normal';
xlabel("c_A")
ylabel("c_C")
colorbar
saveas(gcf,"opt_heat_map_75beds_lognorm.png")
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
m=75;
T = @(k) 1/sqrt(1+k);
X_opt = zeros(100,2);
for j = 1:100
cA0 = randi(25);

X = [cA0,randi(m-cA0)];

i=0;
flag = true;
while flag;
    i =i+1;
    Y = [100,100];
    %%%%%Random walk%%%%%%
    while sum(Y)>m || any(Y<0)
         dirct = mod(i,2); %systematic
         delta_X = [sign(randn)*dirct,sign(randn)*(1-dirct)];
         Y = X(i,:)+delta_X;
    end
    %%%%%end of permutation%%%
    UY = f(Y(1),Y(2));
    UX = f(X(i,1),X(i,2));

    if UY<UX
        X(i+1,:) = Y; %accept 100%
    else
        if rand()<= exp(-(UY-UX)/T(i)) %accept with probability
            X(i+1,:) = Y; 
        else
            X(i+1,:) = X(i,:); %reject
        end
    end
if i>20    
    flag = sum(abs(diff(X(end-10:end,:))),'all')>0;
end
end
X_opt(j,:) = X(end,:);
j
end


function S = f(capA,capC)

% Objective function
% A
%lambda1 = @(t) (-(1/3650).*t.^2 + (1/10).*t); % arrival rate, ward A
mu1 = log(4*sqrt(2));
s2_1 = log(2); 
% Length of stay for A is lognormal dist. 
% Mean and sd of 8 days.

% B
%lambda2 = @(t) ((1/5).*lambda1(t));
mu2 = log(6*sqrt(2));
s2_2 = log(2); 
% Length of stay for B is lognormal dist. 
% Mean and sd of 12 days.

% C
%lambda3 = @(t) (6);
mu3 = log(5*sqrt(2));
s2_3 = log(2);
% Length of stay for C is lognormal dist. 
% Mean and sd of 10 days.
[Rejec, Realloc, ~, ~] = BedUtil([capA,75-capA-capC,capC],[mu1,mu2,mu3],[s2_1, s2_2, s2_3]);
S = sum(Rejec,'all')+sum(Realloc,'all');
end