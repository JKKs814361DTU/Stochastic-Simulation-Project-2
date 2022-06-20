function [Rejec, Realloc, bedocc, no_patients] = BedUtil(cap,mu,sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bedocc = zeros(3,365+1); % Number of beds occupied in each ward
Rejec = zeros(3,365+1); %Number of rejected for ward A and C, B is 0
Realloc = zeros(1,365+1); %Number reloctacted from B to A
no_patients = zeros(1,3); % Needed to determine fraction of rejected pt's.

arr_t = zeros(3,365);
lambda1 = @(t) (-(1/3650)*t^2 + (1/10)*t); % arrival rate, ward A
lambda2 = @(t) ((1/5)*lambda1(t)); %\% arrival rate, ward B
lambda3 = @(t) 6;% arrival rate, ward C
%mean_occ = zeros(3,1);

%Vectors of how long each bed is going to be occupied for. 0 means not
%occupied. E.g. [17,3,0] means bed first bed occupied for 17 more days,
%bed 2 occupied for 3 more days, bed 3 empty
Stay1 = zeros(1,cap(1));  %Beds for ward A
Stay2 =zeros(1,cap(2));%Beds for ward B
Stay3 = zeros(1,cap(3)); %Beds for ward C


%%%Simulation starts%%%%%
for t = 1:365
    pA = round(lambda1(t)); % number of patients arriving day t for A
    pB = round(lambda2(t)); % number of patients arriving day t for A
    pC = round(lambda3(t)); % number of patients arriving day t for A
    
    % counting total number of patient types
    no_patients(1) = no_patients(1) + pA; 
    no_patients(2) = no_patients(2) + pB;
    no_patients(3) = no_patients(3) + pC;
    
    %%% Intensive care patients are first priority (Ward B)
    if (pB <= cap(2)-nnz(Stay2)) %new arrivals to B is less than number of avaliable beds
        bedocc(2,t) = nnz(Stay2) + pB; %add to ward B

    elseif (pB > cap(2)-nnz(Stay2)) %new arrivals are higher than number of avaliable beds
        %if (cap(2)-bedocc(2,t) > 0) %new arrivals are Redundant?
        pA = pA + (pB - (cap(2)-nnz(Stay2))); % redirect  to A
        Realloc(t) = (pB - (cap(2)-nnz(Stay2))); %save number of relocated from B to A
        pB = pB - (cap(2)-nnz(Stay2)); %update new arrivals to B
        bedocc(2,t) = nnz(Stay2) + pB; %add to ward B
        
    end 
    
    %%%Update Ward A
    if (pA <= cap(1)-nnz(Stay1)) %new arrivals to A is less than number of avaliable beds
        bedocc(1,t) = nnz(Stay1) + pA; %add to ward A
    elseif (pA > cap(1)-nnz(Stay1)) %new arrivals are higher than number of avaliable beds
        Rejec(1,t) = (pA - (cap(1)-nnz(Stay1))); % redirect to different hospital
        pA = cap(1)-nnz(Stay1); %arrivals to A that can be accepted
        bedocc(1,t) = nnz(Stay1) + pA; %add to ward A
    end 
    %%%Update Ward C
    if (pC <= cap(3)-nnz(Stay3)) %new arrivals to A is less than number of avaliable beds
        bedocc(3,t) = nnz(Stay3) + pC; %add to ward A
    elseif (pC > cap(3)-nnz(Stay3)) % more arrivals than available beds
        Rejec(3,t) = (pC - (cap(3)-nnz(Stay3))); % redirect to different hospital
        pC = cap(3)-nnz(Stay3); % arrivals to A that can be accepted
        bedocc(3,t) = nnz(Stay3) + pC; %add to ward A
    end 

    %%%%%%%%%%%Departures%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %One day at hospital; Decrease counter of days by 1
    Stay1 = Stay1 - 1;
    Stay2 = Stay2 - 1;
    Stay3 = Stay3 - 1;
    %Correct negative numbers and sort
    Stay1 = sort(max(Stay1,0),'descend');
    Stay2 = sort(max(Stay2,0),'descend');
    Stay3 = sort(max(Stay3,0),'descend');


    %%%Fil up beds
    dA = lognrnd(mu(1),sigma(1),1,pA); %Draw random days of stay for pA new patients
    dB = lognrnd(mu(2),sigma(2),1,pB); %Draw random days of stay for pB new patients 
    dC = lognrnd(mu(3),sigma(3),1,pC); %Draw random days of stay for pC new patients 

    %%Assign day counters to beds
    Stay1(min(find(Stay1==0)):(min(find(Stay1==0))+pA-1)) = dA; 
    Stay2(min(find(Stay2==0)):(min(find(Stay2==0))+pB-1)) = dB;
    Stay3(min(find(Stay3==0)):(min(find(Stay3==0))+pC-1)) = dC;

    bedocc(1:3,t) = [nnz(Stay1);nnz(Stay2);nnz(Stay3)]; %update the number of beds occupied
    
end

end

