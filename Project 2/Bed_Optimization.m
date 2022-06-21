%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utilization of hospital beds during epidemics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
close all

m=75; %%number of beds
iter =1;
%% MATLAB Optimization Toolbox Patternsearch
tic
PSoptions = optimoptions('patternsearch','Display','iter');
optimoptions("patternsearch",'Display','iter')
Objfcn = @(x) f(round(x(1)),round(x(2)),m);
for i =1:iter
cA0 = randi(25);
X0 = [cA0,randi(m-cA0)];
[Xps(i,:),Fps] = patternsearch(Objfcn,X0,[],[],[],[],[1,1],[m,m],PSoptions);
end
toc
figure;

histogram2(Xps(:,1),Xps(:,2),'BinMethod','integers','DisplayStyle','tile','ShowEmptyBins','on')
title("Pattern search optimal solution - "+string(m)+" beds")
xlabel("c_A")
ylabel("c_C")
colorbar
saveas(gcf,"PS_"+string(m)+"_beds.png")
saveas(gcf,"PS_"+string(m)+"_beds")
%% Brute force method

tic
S_10 = zeros(m,m,10);
for i = 1:iter
for capA = 1:m
    for capC =1:m
        if capA+capC <=m
            S_10(capA,capC,i)=f(capA,capC,m);
        else
            S_10(capA,capC,i) =NaN;
        end
    end
end
i
end
toc
%%
S = mean(S_10,3);

% %%
% load Optimiziation_75_beds.mat
% [Cap_A,Cap_C] = find(S==min(S,[],'all'));
% Cap_B=m-Cap_A-Cap_C;
% figure;
% imagesc(S)
% ax = gca;
% ax.YDir = 'normal';
% xlabel("c_A")
% ylabel("c_C")
% colorbar
% saveas(gcf,"opt_heat_map_75beds_lognorm.png")
% figure;
% hold on
% contour(S,'ShowText','on')
% plot(Cap_A,Cap_C,'*r')
% for i = 1:10
% [Cap_A,Cap_C] = find(S_10(:,:,i)==min(S_10(:,:,i),[],'all'));
% plot(Cap_A,Cap_C,'*b')
% end
% 
% xlabel("c_A")
% ylabel("c_C")
% saveas(gcf,"opt_contour_map_75beds_lognorm.png")
% hold off
% 
% save Optimiziation_75_beds.mat


%% Simulated annealing

T = @(k) 1/sqrt(1+k);
X_opt = zeros(iter,2);
tic
for j = 1:iter
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
    UY = f(Y(1),Y(2),m);
    UX = f(X(i,1),X(i,2),m);

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
toc
figure;
title("Simulated annealing optimal solution - " +string(m)+" beds")
histogram2(X_opt(:,1),X_opt(:,2),'BinMethod','integers','DisplayStyle','tile','ShowEmptyBins','on')
xlabel("c_A")
ylabel("c_C")
colorbar
% saveas(gcf,"SA_"+string(m)+"_beds.png")
% saveas(gcf,"SA_"+string(m)+"_beds")
% 
% save("Opt_"+string(m)+"beds.mat")
function S = f(capA,capC,m)

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
if capA + capC >m
    S = inf;
else
[Rejec, Realloc, ~, ~] = BedUtil([capA,m-capA-capC,capC],[mu1,mu2,mu3],[s2_1, s2_2, s2_3]);
S = sum(Rejec,'all')+sum(Realloc,'all');
end
end