clc;clear;
% A
lambda1 = @(t) (-(1/3650).*t.^2 + (1/10).*t); % arrival rate, ward A
mu1 = log(4*sqrt(2));
s2_1 = log(2); 
% Length of stay for A is lognormal dist. 
% Mean and sd of 8 days.

d1 = makedist("Lognormal","mu",mu1,sigma=s2_1);
d2 = makedist("Exponential","mu",mu1);
x=linspace(0,20);
plot(x,pdf(d1,x),x,pdf(d2,x))