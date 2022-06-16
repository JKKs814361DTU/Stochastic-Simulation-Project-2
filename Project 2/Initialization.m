%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utilization of hospital beds during epidemics
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
t0 = 0; t_final = 365;
% Initial guess for distribution of the beds
Cap_A = 10;
Cap_B = 10;
Cap_C = 55;

% A
lambda1 = @(t) (-(1/3650)*t^2 + (1/10)*t); % arrival rate, ward A
mu1 = log(4*sqrt(2));
s2_1 = log(2); 
% Length of stay for A is lognormal dist. 
% Mean and sd of 8 days.

% B
lambda2 = @(t) ((1/5)*lambda1(t));
mu2 = log(6*sqrt(2));
s2_2 = log(2); 
% Length of stay for B is lognormal dist. 
% Mean and sd of 12 days.

% C
lambda3 = 6;
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

[Rejected, Reallocated, mean_occ] = BedUtil([Cap_A, Cap_B, Cap_C],...
    [mu1,mu2,mu3],[s2_1, s2_2, s2_3],[lambda1, lambda2, lambda3]);
