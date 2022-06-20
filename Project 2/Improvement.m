%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Improving the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
rng(19);
% Initial guess for distribution of the beds
% Cap_A = 25;
% Cap_B = 25;
% Cap_C = 25;

% We now need a better guess
Cap_A = 29;
Cap_B = 12;
Cap_C = 29;

mu = [log(4*sqrt(2)) log(6*sqrt(2)) log(5*sqrt(2))];
s = [log(2) log(2) log(2)];

%% Simulation
[Rejected, Reallocated, bedocc, no_patients] = BedUtil([Cap_A, Cap_B, Cap_C],...
    mu,s);

%% 
disp(no_patients)
disp([sum(Rejected(1,:)) sum(Reallocated), sum(Rejected(3,:))])

mnA = sum(Rejected(1,:))/no_patients(1)
mnB = sum(Reallocated(1,:))/no_patients(2)
mnC = sum(Rejected(3,:))/no_patients(3)

%% 
figure();
plot(1:366,bedocc(1,1:end),1:366,bedocc(2,1:end),1:366,bedocc(3,1:end))
legend("A","B","C")
xlabel('days')
title('Beds occupied in each ward, Simulation 2')

figure()
plot(1:366, Rejected(1,1:end), 1:366, Reallocated(1:end), 1:366, Rejected(3,1:end))
title("Rejections and reallocations, Simulation 2")
xlabel('days')
legend("A","B","C")