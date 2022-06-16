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

