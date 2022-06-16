function [Rejec, Realloc, bedocc] = BedUtil(cap,mu,sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bedocc = zeros(3,365+1); % Number of beds occupied in each ward
Rejec = zeros(3,365+1);
Realloc = zeros(3,365+1);

arr_t = zeros(3,365);
lambda1 = @(t) (-(1/3650)*t^2 + (1/10)*t); % arrival rate, ward A
lambda2 = @(t) ((1/5)*lambda1(t));
lambda3 = @(t) 6;
%mean_occ = zeros(3,1);
Stay1 = zeros(1,cap(1));
Stay2 =zeros(1,cap(2));
Stay3 = zeros(1,cap(3));
for t = 1:365
    pA = lambda1(t); % number of patients arriving day t for A
    pB = lambda2(t); % number of patients arriving day t for A
    pC = lambda3(t); % number of patients arriving day t for A
    
    if (pB <= cap(2)-bedocc(2,t)) %new arrivals to B are lower avaliable places
        bedocc(2,t) = bedocc(2,t) + pB; %add to ward B
    elseif (pB > cap(2)-bedocc(2,t)) %new arrivlas are higher than avaliable places
        %if (cap(2)-bedocc(2,t) > 0) %new arrivals are Rundandt?
        pA = pA + (pB - (cap(2)-bedocc(2,t))); % redirect  to A
        pB = cap(2)-bedocc(2,t);
        Reallocated(2,t) = (pB - (cap(2)-bedocc(2,t)));
        bedocc(2,t) = bedocc(2,t) + pB; %add to ward B
        
    end 
    %%%Update Ward A
    if (pA <= cap(1)-bedocc(1,t)) %new arrivals to A are lower avaliable places
        bedocc(1,t) = bedocc(1,t) + pA; %add to ward A
    elseif (pA > cap(1)-bedocc(1,t)) %new arrivlas are higher than avaliable places
        Rejec(t) = pA + (pA - (cap(1)-bedocc(1,t))); % redirect  to reject
        pA = cap(1)-bedocc(1,t); %arrivals to A that can be accepted
        bedocc(1,t) = bedocc(1,t) + pA; %add to ward A
    end 
    %%%Update Ward C
    if (pC <= cap(3)-bedocc(3,t)) %new arrivals to A are lower avaliable places
        bedocc(3,t) = bedocc(3,t) + pC; %add to ward A
    elseif (pC > cap(3)-bedocc(3,t)) %new arrivlas are higher than avaliable places
        Rejec(t) = pC + (pC - (cap(3)-bedocc(3,t))); % redirect  to reject
        pC = cap(3)-bedocc(3,t); %arrivals to A that can be accepted
        bedocc(3,t) = bedocc(3,t) + pC; %add to ward A
    end 

    %%%%%%%%%%%Departures%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dA = lognrnd(mu(1),sigma(1));
    dB = lognrnd(mu(2),sigma(2));
    dC = lognrnd(mu(3),sigma(3));

    bedocc(1:3,t) = bedocc(1:3,t)-[dA;dB;dC];
    bedocc(1:3,t) = max(bedocc(1:3,t),0);

%         if (pB > 0 && cap(1)-bedocc(1,t)>0)
%             if pB < cap(1)-bedocc(1,t)
%                 bedocc(1,t) = bedocc(1,t) + pB;
%                 pB = 0;
%             else if (pB > cap(1) - bedocc(1,t))
%                 bedocc(1,t) = bedocc(1,t) + (cap(1)-bedocc(1,t));
%                 pB = pB - (cap(1)-bedocc(1,t));
%             end
%     if (pA < bedocc(1,t))
    
end

end

