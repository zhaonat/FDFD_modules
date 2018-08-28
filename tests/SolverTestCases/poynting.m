% EE256:Numerical Electromagnetics, Spring 2013-2014
% Homework 4
% Problem 3: Interaction of a plane wave with a scatterer
% Script by Chunliang Zheng

function [Sx, Sy] = poynting(Hz, Ex, Ey)
%% Input Parameters
% Hz, Ex, Ey: 2D arrays of H- and E-field components

%% Output Parameters
% Sx, Sy: 2D array of x- and y-components of Poynting vector
Sx=1/2*real(Ey.*conj(bwdmean_w(Hz,'x')));
Sy=-1/2*real(Ex.*conj(bwdmean_w(Hz,'y')));


%% Grading breakdown: (total 1 point)

% You get 1 point if you calculate the Poynting vector correctly