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

for i = 1:size(Cap,1)
    [Rejected, Reallocated, bedocc, no_patients] = BedUtil([Cap(i,1), Cap(i,2), Cap(i,3)],...
    [mu1,mu2,mu3],[s2_1, s2_2, s2_3]);
    
    figure()
    subplot(1,2,1)
        plot(1:366, Rejected(1,1:end), 1:366, Reallocated(1:end), 1:366, Rejected(3,1:end))
        title("Redirected and relocated patients");
        xlabel('days')
        legend("A","B","C")
    subplot(1,2,2)    
        plot(1:366,bedocc(1,1:end),1:366,bedocc(2,1:end),1:366,bedocc(3,1:end))
        legend("A","B","C")
        xlabel('days')
        title("Beds occupied in each ward");
    sgtitle("Total capacity: " + sum(Cap(i,:)));
    
    disp("Capacity: " + sum(Cap(i,:)))
    disp(no_patients)
    disp([sum(Rejected(1,:)) sum(Reallocated), sum(Rejected(3,:))])
    disp([sum(Rejected(1,:))/no_patients(1) sum(Reallocated)/no_patients(2), sum(Rejected(3,:))/no_patients(3)])
    disp('Mean fraction of beds occupied in each ward');
    disp([mean(bedocc(1,:)/Cap(i,1)), mean(bedocc(2,:)/Cap(i,2)), mean(bedocc(3,:)/Cap(i,3)) ])

end
