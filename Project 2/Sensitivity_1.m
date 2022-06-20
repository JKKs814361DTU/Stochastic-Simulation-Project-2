%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utilization of hospital beds during epidemics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
close all

%% Initialization
rng(19);
% Initial guess for distribution of the beds
Cap = [20 20 20; 27 26 27; 34 33 33];

mu1 = log(4*sqrt(2));
s2_1 = log(2); 
mu2 = log(6*sqrt(2));
s2_2 = log(2); 
mu3 = log(5*sqrt(2));
s2_3 = log(2);

for i = 1:size(Cap,2)
    [Rejected, Reallocated, bedocc, no_patients] = BedUtil([Cap(i,1), Cap(i,2), Cap(i,3)],...
    [mu1,mu2,mu3],[s2_1, s2_2, s2_3]);

    disp(no_patients)
    disp([sum(Rejected(1,:)) sum(Reallocated), sum(Rejected(3,:))])

    mnA = sum(Rejected(1,:))/no_patients(1)
    mnB = sum(Reallocated(1,:))/no_patients(2)
    mnC = sum(Rejected(3,:))/no_patients(3)

end
    


%% Simulation

% Primary performance measure:
% Probability that all beds are occupied on arrival for each of the three
% patient types as well as the mean number of patients, that are relocated
% due to shortage of beds.
% Also estimate the mean fraction of beds that are utilized in each ward.



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

% Mean fraction of beds occupied in each ward
figure();
plot(1:366,bedocc(1,1:end)/Cap_A,1:366,bedocc(2,1:end)/Cap_B,1:366,bedocc(3,1:end)/Cap_C)
legend("A","B","C")
xlabel('days')
title('Fraction of beds occupied in each ward')

%[bedocc(1,:)/Cap_A, bedocc(2,:)/Cap_B, bedocc(3,:)/Cap_C];
disp('Mean fraction of beds occupied in each ward')
disp([mean(bedocc(1,:)/Cap_A), mean(bedocc(2,:)/Cap_B), mean(bedocc(3,:)/Cap_C) ])


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