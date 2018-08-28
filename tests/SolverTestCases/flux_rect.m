% EE256:Numerical Electromagnetics, Spring 2013-2014
% Homework 4
% Problem 3: Interaction of a plane wave with a scatterer
% Script by Chunliang Zheng

function P = flux_rect(Sx, Sy, ind_x, ind_y, dL)
%% Input Parameters
% Sx, Sy: 2D arrays of Poynting vectors in x- and y-directions
% ind_x, ind_y: x- and y-indices of rectangular region
% dL: [dx dy]

%% Output Parameters
% P: net power flux going out of rectangular region
Pxn=sum(Sx(ind_x(1),ind_y(1):ind_y(2)))*dL(2);
Pyn=sum(Sy(ind_x(1):ind_x(2),ind_y(1)))*dL(1);

Pxp=sum(Sx(ind_x(2)+1,ind_y(1):ind_y(2)))*dL(2);
Pyp=sum(Sy(ind_x(1):ind_x(2),ind_y(2)+1))*dL(1);

P=Pxp-Pxn+Pyp-Pyn;


%% Grading breakdown: (total 1 point)

% You get 1 point if you calculate the net power flux correctly