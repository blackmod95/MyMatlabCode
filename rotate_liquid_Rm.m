clc;
clear;

alpha = 4;
k = 9;
H0 = 20;
g = 9.8;
mu = 2;
ro = 3;
b0 = 4;

sig1 = 0.1:0.1:10;
sig2 = 0;
sigma = sig1 + j*sig2;

sigma4 = sigma.^4
sigma3 = sigma.^3
sigma2 = sigma.^2

A = sigma4 + (g*H0*k^2 - alpha^2 - 2*b0*k^2/(mu*ro)).*sigma2 - (g*H0*b0^2*k^4)/(mu*ro) + (b0*k)^4/(mu*ro)^2;
B = j*2*k^2.*sigma3 + j*(g*H0*k^4-2*alpha^2*k^2).*sigma + j*(2*b0^2*k^4.*sigma)/(mu*ro);
C = -2*k^4.*sigma2 + alpha^2*k^4;

R_m = (-B - sqrt(B.^2 - 4*A.*C))./(2*A);

[Rm, idxs] = sort(R_m);

plot(Rm, real(sigma(idxs)));
