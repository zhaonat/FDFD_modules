%% derivative operators with dirichlet boundary conditions
close all;
clear; 

% Nx = 100; 
% Ny = 100; 
% Nz = 10;
% N = [Nx, Ny, Nz]; dx = 0.01; dy = dx; dz = dx;
% dL = [dx dy, dz];

Nx = 100; 
Ny = 100; 
N = [Nx, Ny]; dx = 0.01; dy = dx; dz = dx;
dL = [dx dy];
%% Solver Code Begins Here
w = 'x'; s = 'b';
%w = 'x', 'y', or 'z'
%s = 'b' or 'w'
% dL: [dx dy dz] for 3D; [dx dy] for 2D
% N: [Nx Ny Nz] for 3D; [Nx Ny] for 2D

dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
M = prod(N);  % total number of cells in domain
%take advantage of linear indexing
dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
M = prod(N);  % total number of cells in domain
%take advantage of linear indexing

ind_cur = 1:M;  % indices of current points
ind_cur = ind_cur(:);

ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
ind_adj_x = reshape(ind_adj_x, N);
ind_adj_x0 = ind_adj_x;

%% shift the row indices
Dws = -sign*(5/2)*speye(M); %M fully determines the matricial size;
for w = ['x', 'y']
    ind_adj_x = circshift(ind_adj_x0,  -sign * ('xyz' == w));
    ind_adj_2x = circshift(ind_adj_x0, -2*sign * ('xyz' == w));
    ind_adj_rx = circshift(ind_adj_x0,  sign * ('xyz' == w));
    ind_adj_r2x = circshift(ind_adj_x0, 2*sign * ('xyz' == w));

    ind_adj_xl = ind_adj_x(:);
    ind_adj_2xl = ind_adj_2x(:);
    ind_adj_rxl = ind_adj_rx(:);
    ind_adj_r2xl = ind_adj_r2x(:);

    %% conver the offdiagonals into a linear index
    linear_ind_x = sub2ind([M M], ind_cur, ind_adj_xl); %points = rowsub...ind_adj = colsub
    linear_ind_2x = sub2ind([M M], ind_cur, ind_adj_2xl); %points = rowsub...ind_adj = colsub
    linear_ind_rx = sub2ind([M M], ind_cur, ind_adj_rxl); %points = rowsub...ind_adj = colsub
    linear_ind_r2x = sub2ind([M M], ind_cur, ind_adj_r2xl); %points = rowsub...ind_adj = colsub
    
    Dws(linear_ind_x) = (4/3)*sign;
    Dws(linear_ind_2x) = (-1/12)*sign;
    Dws(linear_ind_rx) = (4/3)*sign;
    Dws(linear_ind_r2x) = (-1/12)*sign;

end
Dws = (1/dw^2)*Dws;
omega = 155;
constant = 1e3
Dxf5 = createDws_dense('x','f', dL,N);
Dxb5 = createDws_dense('x','b', dL,N);
Dyf5 = createDws_dense('y','f', dL,N);
Dyb5 = createDws_dense('y','b', dL,N);
Ax5 = Dxf5*Dxb5; Ay5 = Dyf5*Dyb5;

%% apparently, the frequency mapping is different between the 9 point and the
%% 5 point stencil, which is odd in a physical sense since they are the 
%% same problem, numerically, it doesn't because A5 and A are different
A5 = Ax5+Ay5 + constant*speye(N.^2);
A = Dws + 23.99*constant*speye(N.^2);
% constant has to be tuned

b = zeros(length(Dws),1);
b(4750) = 1i*omega

x1 = A\b;
x2 = A5\b;
E = reshape(x1, N(1), N(2));
E5 = reshape(x2, N(1), N(2));
subplot(1,2,1)
visabs(E, [-1,1], [-1,1]);
subplot(1,2,2)
visabs(E5, [-1,1], [-1,1]);
