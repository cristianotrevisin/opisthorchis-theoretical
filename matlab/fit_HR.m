clc
close all
clear all

x = 1:10000; x = x';
y = 20*x.^(-0.82).*x;

ft = fittype( 'hyperbolic_relation( x, rho, omega )')
f = fit(x, y, ft)
figure
plot(f,x,y)