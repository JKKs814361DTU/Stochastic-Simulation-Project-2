%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation using the optimal distribution of 75 beds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc
close all

%% Initialization
rng(19);
Cap = [26 14 35];
Total_Beds = sum(Cap);

mu1 = log(4*sqrt(2));
s2_1 = log(2); 
mu2 = log(6*sqrt(2));
s2_2 = log(2); 
mu3 = log(5*sqrt(2));
s2_3 = log(2);

[Rejected, Reallocated, bedocc, no_patients] = BedUtil([Cap(1), Cap(2), Cap(3)],...
[mu1,mu2,mu3],[s2_1, s2_2, s2_3]);

figure()
subplot(1,2,1)
    plot(1:366, Rejected(1,1:end), 1:366, Reallocated(1:end), 1:366, Rejected(3,1:end))
    title("Redirections and relocations");
    xlabel('days')
    legend("A","B","C")
subplot(1,2,2)    
    plot(1:366,bedocc(1,1:end),1:366,bedocc(2,1:end),1:366,bedocc(3,1:end))
    legend("A","B","C")
    xlabel('days')
    title("Beds occupied in each ward");
sgtitle("Optimal distribution of " + sum(Cap(:)) + " beds");

disp("Capacity: " + sum(Cap(:)))
disp(no_patients)
disp([sum(Rejected(1,:)) sum(Reallocated), sum(Rejected(3,:))])
disp([sum(Rejected(1,:))/no_patients(1) sum(Reallocated)/no_patients(2), sum(Rejected(3,:))/no_patients(3)])
disp('Mean fraction of beds occupied in each ward');
disp([mean(bedocc(1,:)/Cap(1)), mean(bedocc(2,:)/Cap(2)), mean(bedocc(3,:)/Cap(3)) ])

