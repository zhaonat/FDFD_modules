time1=fix(clock);                   % start time, [sec]
c0 = 299792458;                   % speed of light in vacuum, [m/s]
e0=8.854187817e-12;         % vacuum permittivity 
mu0=pi*4e-7;         

epsilon_0 = 1;
epsilon_1 = 12;

kx = 0; % normal indidence

lambda_scan = linspace(1, 2, 1000)*1e-6;
omega_scan = 2*pi*c0./lambda_scan;

kz0 = sqrt(epsilon_0*omega_scan.^2/c0^2 - kx^2);
kz1 = sqrt(epsilon_1*omega_scan.^2/c0^2 - kx^2);
h = 0.76*1e-6;

c1 = (kz0.^2-kz1.^2)./(2*kz1.*kz0);
c2 = (kz0.^2+kz1.^2)./(2*kz1.*kz0);
r_d = 1i*c1.*sin(kz1*h)./(cos(kz1*h)-1i*c2.*sin(kz1*h));

t_d = 1./((cos(kz1*h)-1i*c2.*sin(kz1*h)));

figure();
plot(lambda_scan, abs(r_d).^2);
hold on;
plot(lambda_scan, abs(t_d).^2);
hold on;
plot(lambda_scan, abs(r_d).^2+abs(t_d).^2)