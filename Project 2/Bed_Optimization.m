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
xlabel("c_A")
ylabel("c_C")
saveas(gcf,"opt_contour_map_75beds_lognorm.png")
hold off

save Optimiziation_75_beds.mat
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