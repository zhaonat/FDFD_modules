%% sfactor test cases
%% CREATE SFACTOR ALREADY ANTICIPATES THAT YOU HAVE CREATED A CONSISTENT DOMAIN PML
%%===============0 length PML========================00

wrange = [-6*10^-6 6*10^-6]
s = 'f'
c = 3*10^8; c0 = c;
eps0 = 8.85*10^-12;
mu0 = 4*pi*10^-7;
Nw = 100
Nw_pml = 10;
Nw = Nw+Nw_pml; %this is always true in the actual code;

L0 = 10^-6; wvlen = 3*L0;
omega = 2*pi*c0/(wvlen);
lnR = -16;

sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml);
figure;
plot(abs(sfactor_array));

%%===================Excessive Length PML======================

Nw_pml = 10;
Nw = Nw+Nw_pml;
sfactor_array = create_sfactor_mine(wrange, s, omega, eps0, mu0, Nw, Nw_pml)
figure;
plot(abs(sfactor_array));

%%================Compare to hw solutoin==================%
sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml)
figure;
plot(abs(sfactor_array));

