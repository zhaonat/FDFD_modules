%based on Kumar et al. Hyperbolic materials made


%h-BN
function [epsilonv,epsilonh]=hBNNEW(lambda) % lambda is in micron

c0=2.99792458e8;
w=2*pi*c0./(lambda*1e-6);

wn=1./lambda*1e4; % in cm-1


e_infh=4.87; e_infv=2.95;

rh=5; rv=4;

wnLOA2u=830; wnTOA2u=780; %this is a vertical mode
wnLOE1u=1610; wnTOE1u=1370; %this is a horizontal mode

epsilonv=e_infv+e_infv*(        (wnLOA2u^2-wnTOA2u^2)./(wnTOA2u^2-1i*rv*wn-wn.^2)       ); %+      (wnLOE1u^2-wnTOE1u^2)./(wnTOE1u^2-1i*rv*wn-wn.^2)

epsilonh=e_infh+e_infh*(         (wnLOE1u^2-wnTOE1u^2)./(wnTOE1u^2-1i*rh*wn-wn.^2)    ); %  (wnLOA2u^2-wnTOA2u^2)./(wnTOA2u^2-1i*rh*wn-wn.^2)    +  

end