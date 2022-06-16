function [Rejected, Reallocated, mean_occ] = BedUtil(cap,mu,sigma,lambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bedoc = zeros(3,365+1); % Number of beds occupied in each ward
Rejec = zeros(3,365+1);
Realloc = zeros(3,365+1);

arr_t = zeros(3,365

mean_occ = zeros(3,1);

for t = 0:365
    pA = lambda(1)(t); % number of patients arriving day t for A
    pB = lambda(2)(t); % number of patients arriving day t for A
    pC = lambda(3)(t); % number of patients arriving day t for A
    
    if (pB < cap(2)-bedocc(2,t))
        bedocc(2,t) = bedocc(2,t) + pB;
    else if (pB > cap(2)-bedocc(2,t))
        if (cap(2)-bedocc(2,t) > 0)
            bedocc(2,t) = bedocc(2,t) + cap(2)-bedocc(2,t);
            pB = pB - cap(2)-bedocc(2,t);
        end 
        if (pB > 0 && cap(1)-bedocc(1,t)>0)
            if pB < cap(1)-bedocc(1,t)
                bedocc(1,t) = bedocc(1,t) + pB;
                pB = 0;
            else if (pB > cap(1) - bedocc(1,t))
                bedocc(1,t) = bedocc(1,t) + (cap(1)-bedocc(1,t));
                pB = pB - (cap(1)-bedocc(1,t));
            end
    if (pA < bedocc(1,t))
    
end

end

