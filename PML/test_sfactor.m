%% sfactor test cases
%% CREATE SFACTOR ALREADY ANTICIPATES THAT YOU HAVE CREATED A CONSISTENT DOMAIN PML
%%===============0 length PML========================00

wrange = [-6 6]
s = 'f'
eps0 = 8.85*10^-12*L0;
mu0 = 4*pi*10^-7*L0;
c = 1/sqrt(mu0*eps0); c0 = c;

Nw = 100;
Nw_pml = 10;
Nw = Nw+Nw_pml; %this is always true in the actual code;

L0 = 10^-6; wvlen = 3;
omega = 2*pi*c0/(wvlen);
lnR = -16;

sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml);
figure;
plot(abs(sfactor_array));

%%===================Excessive Length PML======================

sfactor_array_mine = create_sfactor_mine(wrange, s, omega, eps0, mu0, Nw, Nw_pml);
figure;
plot(abs(sfactor_array_mine));

%%================Compare to hw solutoin==================%
figure;
plot(abs(sfactor_array));
hold on;
plot(abs(sfactor_array_mine))

