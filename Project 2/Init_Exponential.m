%% Initialization using an exponential distribution for the LOS 
% pretty much the same initialization, but with minor changes for the run
% of LOS being exponential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc 
close all
%% Introduction

% A: Regular care ward, 
% B: Intensive care ward
% C: Inpatient ward, originally contained 75 beds, some of which are now
% moved to A and B.

% t0 = 0,    t_final = 365.

%% Initialization
rng(19);
% Initial guess for distribution of the beds
Cap_A = 25;
Cap_B = 25;
Cap_C = 25;

% better guess for initial distribution of beds? 

% A
lambda1 = @(t) (-(1/3650).*t.^2 + (1/10).*t); % arrival rate, ward A
mu1 = log(4*sqrt(2));
s2_1 = log(2); 
% Length of stay for A is lognormal dist. 
% Mean and sd of 8 days.

% B
lambda2 = @(t) ((1/5).*lambda1(t));
mu2 = log(6*sqrt(2));
s2_2 = log(2); 
% Length of stay for B is lognormal dist. 
% Mean and sd of 12 days.

% C
lambda3 = @(t) (6);
mu3 = log(5*sqrt(2));
s2_3 = log(2);
% Length of stay for C is lognormal dist. 
% Mean and sd of 10 days.


%% Simulation

% Primary performance measure:
% Probability that all beds are occupied on arrival for each of the three
% patient types as well as the mean number of patients, that are relocated
% due to shortage of beds.
% Also estimate the mean fraction of beds that are utilized in each ward.

[Rejected, Reallocated, bedocc, no_patients] = bedUtil_Exponential([Cap_A, Cap_B, Cap_C],...
    [mu1,mu2,mu3],[s2_1, s2_2, s2_3]);

%% Determining fractions
% Estimate the probability that all beds are occupied on arrival for each
% of the three patient types as well as the mean number of patients that
% are relocated due to shortage of beds. Estimate the latter for each
% individual patient type and as a sum over all types. Furthermore,
% estimate the mean fraction of beds that are utilized (occupied) in each
% ward.

disp(no_patients)
disp([sum(Rejected(1,:)) sum(Reallocated), sum(Rejected(3,:))])

mnA = sum(Rejected(1,:))/no_patients(1)
mnB = sum(Reallocated(1,:))/no_patients(2)
mnC = sum(Rejected(3,:))/no_patients(3)


%%
% for i =1:366
% 
% plot(1:i,bedocc(1,1:i),1:i,bedocc(2,1:i))
% legend("A","B","C")
% %pause(0.25)
% end
figure();
plot(1:366,bedocc(1,1:end),1:366,bedocc(2,1:end),1:366,bedocc(3,1:end))
legend("A","B","C")
xlabel('days')
title('Beds occupied in each ward')
%%
figure();
plot(Rejected')
title("Reject")
legend("A","B","C")


figure()
plot(1:366, Rejected(1,1:end), 1:366, Reallocated(1:end), 1:366, Rejected(3,1:end))
title("Rejections and reallocations")
xlabel('days')
legend("A","B","C")
%%
figure;
plot(Reallocated')
%legend("A","B","C")
title("reloc")

%%
figure();
plot(1:366,lambda1(1:366))
hold on
plot(1:366,lambda2(1:366))
yline(lambda3(1),'g')
title('Arrival rates')
legend("A","B","C")
xlabel('days')
ylabel('Number of patients')
xlim([0 366]);
%%
% for i =1:366
% 
% plot(1:i,bedocc(1,1:i),1:i,bedocc(2,1:i))
% legend("A","B","C")
% %pause(0.25)
% end
figure();
plot(1:366,bedocc(1,1:end),1:366,bedocc(2,1:end),1:366,bedocc(3,1:end))
legend("A","B","C")
xlabel('days')
title('Beds occupied in each ward')

% Mean fraction of beds occupied in each ward
figure();
plot(1:366,bedocc(1,1:end)/Cap_A,1:366,bedocc(2,1:end)/Cap_B,1:366,bedocc(3,1:end)/Cap_C)
legend("A","B","C")
xlabel('days')
title('Fraction of beds occupied in each ward')

%[bedocc(1,:)/Cap_A, bedocc(2,:)/Cap_B, bedocc(3,:)/Cap_C];
disp('Mean fraction of beds occupied in each ward')
disp([mean(bedocc(1,:)/Cap_A), mean(bedocc(2,:)/Cap_B), mean(bedocc(3,:)/Cap_C) ])