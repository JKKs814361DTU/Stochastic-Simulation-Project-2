%% Sanity check
clear, clc
close all
%%
rng(19);
mu = [log(4*sqrt(2)) log(6*sqrt(2)) log(5*sqrt(2))];
s = [log(2) log(2) log(2)];
cap = [25 25 25];
cap = [32 12 31];
bedocc = zeros(3,365+1); % Number of beds occupied in each ward
Rejec = zeros(3,365+1); %Number of rejected for ward A and C, B is 0
Realloc = zeros(1,365+1); %Number reloctacted from B to A
no_patients = zeros(1,3); % Needed to determine fraction of rejected pt's.

arr_t = zeros(3,365);
lambda1 = @(t) (-(1/3650)*t^2 + (1/10)*t); % arrival rate, ward A
lambda2 = @(t) ((1/5)*lambda1(t)); %\% arrival rate, ward B
lambda3 = @(t) 6;% arrival rate, ward C

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
        pB = pB - Realloc(t); %(cap(2)-nnz(Stay2)); %update new arrivals to B
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
    dA = 8; % mean LOS
    dB = 12; % mean LOS
    dC = 10; % mean LOS

    %%Assign day counters to beds
    Stay1(min(find(Stay1==0)):(min(find(Stay1==0))+pA-1)) = dA; 
    Stay2(min(find(Stay2==0)):(min(find(Stay2==0))+pB-1)) = dB;
    Stay3(min(find(Stay3==0)):(min(find(Stay3==0))+pC-1)) = dC;

    bedocc(1:3,t) = [nnz(Stay1);nnz(Stay2);nnz(Stay3)]; %update the number of beds occupied
    
end

%%
disp('Number of patients of each type:')
disp(no_patients)
disp('Number of redirected or reallocated patients:')
disp([sum(Rejec(1,:)) sum(Realloc), sum(Rejec(3,:))])

mnA = sum(Rejec(1,:))/no_patients(1);
mnB = sum(Realloc(1,:))/no_patients(2);
mnC = sum(Rejec(3,:))/no_patients(3);
disp('Fraction of patients that are redirected:')
disp([mnA, mnB, mnC])

disp('Mean fraction of beds occupied in each ward')
disp([mean(bedocc(1,:)/cap(1)), mean(bedocc(2,:)/cap(2)), mean(bedocc(3,:)/cap(3)) ])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation for different distributions of beds
rng(19);

mu = [log(4*sqrt(2)) log(6*sqrt(2)) log(5*sqrt(2))];
s = [log(2) log(2) log(2)];

capacities = zeros(3,75);
frac_blocked = zeros(75,75,3);
mean_occ = zeros(75,75,3);


%%
x = 0;
X1 = []; Y1 = []; Z1 = [];
X2 = []; Y2 = []; Z2 = [];
X3 = []; Y3 = []; Z3 = [];
X4 = []; Y4 = []; Z4 = [];
X5 = []; Y5 = []; Z5 = [];
X6 = []; Y6 = []; Z6 = [];


capB = 0;
for capA = 74:-1:1
    for capC = 1:(75-capA)
        capB = 75-capA - capB-1;
        if (capA + capB + capC == 75)
            x = x+1;
            cap = [capA capB capC];
            capacities(1,capA) = capA;
            capacities(2,capA) = capB;
            capacities(3,capA) = capC;

            bedocc = zeros(3,365+1); % Number of beds occupied in each ward
            Rejec = zeros(3,365+1); %Number of rejected for ward A and C, B is 0
            Realloc = zeros(1,365+1); %Number reloctacted from B to A
            no_patients = zeros(1,3); % Needed to determine fraction of rejected pt's.

            arr_t = zeros(3,365);
            lambda1 = @(t) (-(1/3650)*t^2 + (1/10)*t); % arrival rate, ward A
            lambda2 = @(t) ((1/5)*lambda1(t)); %\% arrival rate, ward B
            lambda3 = @(t) 6;% arrival rate, ward C

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
                    pB = pB - Realloc(t); %(cap(2)-nnz(Stay2)); %update new arrivals to B
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
                dA = lognrnd(mu(1),s(1),1,pA); %Draw random days of stay for pA new patients
                dB = lognrnd(mu(2),s(2),1,pB); %Draw random days of stay for pB new patients 
                dC = lognrnd(mu(3),s(3),1,pC); %Draw random days of stay for pC new patients 

                %%Assign day counters to beds
                Stay1(min(find(Stay1==0)):(min(find(Stay1==0))+pA-1)) = dA; 
                Stay2(min(find(Stay2==0)):(min(find(Stay2==0))+pB-1)) = dB;
                Stay3(min(find(Stay3==0)):(min(find(Stay3==0))+pC-1)) = dC;

                bedocc(1:3,t) = [nnz(Stay1);nnz(Stay2);nnz(Stay3)]; %update the number of beds occupied

            end
            frac_blocked(capA, capC,1) = sum(Rejec(1,:))/no_patients(1);
            frac_blocked(capA, capC,2) = sum(Realloc(:))/no_patients(2);
            frac_blocked(capA, capC,3) = sum(Rejec(3,:))/no_patients(3);

            X1(x) = capA; Y1(x) = capC; Z1(x) = sum(Rejec(1,:))/no_patients(1);
            X2(x) = capA; Y2(x) = capC; Z2(x) = sum(Realloc)/no_patients(2);
            X3(x) = capA; Y3(x) = capC; Z3(x) = sum(Rejec(3,:))/no_patients(3);

            mean_occ(capA, capC, 1) = mean(bedocc(1,:))/cap(1);
            mean_occ(capA, capC, 2) = mean(bedocc(2,:))/cap(2);
            mean_occ(capA, capC, 3) = mean(bedocc(3,:))/cap(3);

            X4(x) = capA; Y4(x) = capC; Z4(x) = mean(bedocc(1,:))/cap(1);
            X5(x) = capA; Y5(x) = capC; Z5(x) = mean(bedocc(2,:))/cap(2);
            X6(x) = capA; Y6(x) = capC; Z6(x) = mean(bedocc(3,:))/cap(3);
        end
    end
end
%%
subplot(3,2,1);
    plot3(X1, Y1, Z1);
    title('Fraction of rejected patients of type A')
subplot(3,2,2); 
    plot3(X2, Y2, Z2);
    title('Fraction of rejected patients of type B')
subplot(3,2,3); 
    plot3(X3, Y3, Z3);
    title('Fraction of rejected patients of type C')
subplot(3,2,4); 
    plot3(X4, Y4, Z4);
    title('Mean fraction of beds occupied in ward A')
subplot(3,2,5); 
    plot3(X5, Y5, Z5);
    title('Mean fraction of beds occupied in ward B')
subplot(3,2,6); 
    plot3(X6, Y6, Z6);
    title('Mean fraction of beds occupied in ward C')

%%
figure()
    plot3(X1, Y1, Z1);
    title('Fraction of rejected patients of type A')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('frac')
    grid on
figure()
    plot3(X2, Y2, Z2);
    title('Fraction of rejected patients of type B')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('frac')
    grid on
figure()
    plot3(X3, Y3, Z3);
    title('Fraction of rejected patients of type C')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('frac')
    grid on
figure()
    plot3(X4, Y4, Z4);
    title('Mean fraction of beds occupied in ward A')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('mean occ')
    grid on
figure()
    plot3(X5, Y5, Z5);
    title('Mean fraction of beds occupied in ward B')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('mean occ')
    grid on

figure()
    plot3(X6, Y6, Z6);
    title('Mean fraction of beds occupied in ward C')
    xlabel('Cap A'); ylabel('Cap B'); zlabel('mean occ')
    grid on

%%
figure()
hold on
subplot(2,1,1)
    plot3(X1, Y1, Z1);
    subplot(2,1,2);
    %plot3(X2, Y2, Z2);
    %plot3(X3, Y3, Z3);
    plot3(X4, Y4, Z4);
    %plot3(X5, Y5, Z5);
   %plot3(X6, Y6, Z6);
    %legend('Ward A, fraction blocked',...  %'Ward B, fraction blocked', 'Ward C, fraction blocked'
      %  'Ward A, mean occupied frac') % ,'Ward B, mean occupied frac', 'Ward C, mean occupied frac'

%%
plot([1:75],frac_blocked(:,1,1))
hold on
plot([1:75],frac_blocked(:,2,1))
    
%%
hold off
figure()
plot(capacities')
legend('A', 'B', 'C')
sum(capacities,1)