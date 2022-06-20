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

mu = [log(4*sqrt(2)) log(6*sqrt(2)) log(5*sqrt(2))]';
s = [log(2) log(2) log(2)]';

%% Simulation
[Rejected, Reallocated, bedocc, no_patients] = BedUtil([Cap_A, Cap_B, Cap_C],...
    mu,s);