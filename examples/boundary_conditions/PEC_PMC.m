%% tests for the PEC and PMC boundary condition (aka dirichlet for maxwell?)
clear all

% first, we will test the derivative operators (code named Dirichlet)
N = [10,10];
dL = [1,1];
Dxf = createDws_dirichlet('x', 'f', dL,N);
Dxf_2 = createDws_dirichlet('x', 'f', dL,N);
M = prod(N);

% mask creation;
xn = 1:N(1);
yn = 1:N(2);
[Xn,Yn] = meshgrid(xn,yn);
Xn = Xn.'; Yn = Yn.';
mask = ones(N);
mask(Xn == 1) = 0;
mask(Xn == N(1)) =0;

% create an operator mask which simply zeros out elements of the 
% electric field?
figure();
imagesc(mask);
PEC_mask = spdiags(mask(:), 0, M,M);

% if the mask idea is correct, we should be able to reproduce it with the
% dirichlet difference operators...

% create  PEC
Dxfd = createDws_dirichlet_2D('x', 'f', dL,N);
Dxbd = createDws_dirichlet_2D('x', 'b', dL,N);
Dyfd = createDws_dirichlet_2D('y', 'f', dL,N);
Dybd = createDws_dirichlet_2D('y', 'b', dL,N);

del2x = (Dxbd*Dxfd); %% Should the last on-diagonal term be a -2 or a -1?
del2y = Dybd*Dyfd; %% Should the last on-diagonal term be a -2 or a -1?


figure(); spy(Dxbd*Dxfd);
Dxf = createDws('x', 'f', dL,N);
Dxb = createDws('x', 'b', dL,N);
Dyf = createDws('y', 'f', dL,N);
Dyb = createDws('y', 'b', dL,N);

%apparently not...
A = Dxbd*Dxfd+Dyb*Dyf;
Apbc = PEC_mask*Dxb*Dxf+Dyb*Dyf*PEC_mask;

%% the real challenge is actually to know write out maxwell's equations,
% since E and H field follow different differencing operators...

