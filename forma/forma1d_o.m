function [k, D] = forma1d(x, Dt, gamma)
% FORMA1D   1D implementation of FORMA
%
% [k, D] = FORMA1D(x, Dt, gamma) calculates the values of the trap
%           stiffness k and the particle diffusion coefficient D given a
%           sequence of particle positions x sampled with sampling
%           frequency 1/Dt and the friction coefficient gamma of the
%           particle.
%
%   This code is provided with the article:
%
%   High-Performance Reconstruction of Microscopic Force Fields from
%   Brownian Trajectories
%   Laura Perez Garcia, Jaime DonLucas, Giorgio Volpe, Alejandro V. Arzola
%   & Giovanni Volpe 
%   2018

dx = diff(x); 
x = x(1:end-1);
f = gamma*dx/Dt;

k = -sum(x.*f)/sum(x.^2);

D = mean(Dt/(2*gamma^2)*(f+k*x).^2);