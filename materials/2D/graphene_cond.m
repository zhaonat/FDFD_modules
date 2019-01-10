function[sigma_D, sigma_I, sigma] = graphene_cond(w,T,mu,tau)
%tau [unit s], w [unit rad/s], mu[ev], T[K]
%example imput: 1e14,1e-3,0,1e-13

h = 1.054571596e-34;           %Planck's constant over 2*pi
c0 = 2.99792458e+8;            %speed of light in vacuum
kb = 1.3806503e-23;
qe = 1.602176462e-19;          %elementary electric charge %

if T<1
sigma_D=1i/(w+1i/tau)*2*qe^2/pi/h^2*mu*qe/2;
else
sigma_D=1i/(w+1i/tau)*2*qe^2*kb*T/pi/h^2*log(2*cosh(mu*qe/2/kb/T));
%% intra-band term
end

eta=[1e-4:3e-5:1]*2e-17;
for ii=1:length(eta)
    temp(ii)=1i*4*h*w/pi*(G(eta(ii),T,mu)-G(h*w/2,T,mu))/(h^2*w^2-4*eta(ii)^2);
end
sigma_inte=(temp(1)+4*sum(temp(2:2:(end-1)))+2*sum(temp(3:2:(end-2)))+temp(end))*(eta(3)-eta(2))/3;
sigma_I=qe^2/4/h*(G(h*w/2,T,mu)+sigma_inte);
sigma=sigma_D+sigma_I;
end

function Geta=G(eta,T,mu)
kb = 1.3806503e-23;
qe = 1.602176462e-19; 
if(eta/kb/T>700&&eta/kb/T-mu*qe/kb/T>4)
    Geta=1;
elseif(mu*qe/kb/T>700&&mu*qe/kb/T-eta/kb/T>4)
    Geta=0;
else
Geta=sinh(eta/kb/T)/(cosh(mu*qe/kb/T)+cosh(eta/kb/T));
end

end

