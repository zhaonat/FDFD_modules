function Farfield = NTFF_Transform_2D(E_near, H_near, norm_surf, r_near, ...
    dl, r_hat, wavelength, polarization)
% function Farfield = NTFF_Transform_2D(E, H, norm_surf, r_near, dl, r_hat, wavelength, polarization)
% If s-polared, E = Ez [1*n], 
%H = [Hx; Hy] [2*n], if p-polarized, E = {Ex,
% Ey}, H = Hz; norm_surf [2*n], r_near [rx; ry], r_hat is observe direction 
% [2*1]; 
%dl is a scalar... but what should it be?
%H is the normalized H (Z0H). Ignore phasor exp(-jkr)/sqrt(r).

k = 2*pi/wavelength;
factor = exp(1i*pi/4)*sqrt(k/8/pi);

nx = norm_surf(1,:);
ny = norm_surf(2,:);
rx = r_hat(1, :);
ry = r_hat(2, :);
rprime_x = r_near(1,:);
rprime_y = r_near(2,:);

if polarization == 's'
    Ez = E_near;
    Hx = H_near(1,:);
    Hy = H_near(2,:);
    integrand = (- (nx.*Hy - ny.*Hx) + Ez .* (nx.*rx + ny.*ry)) .* exp(1i*k*(rx.*rprime_x + ry.*rprime_y));
    Ez_far = factor * sum(integrand(:)) * dl;
    Farfield = Ez_far;
elseif polarization == 'p'
    Ex = E_near(1,:);
    Ey = E_near(2,:);
    Hz = H_near;
    integrand = ((nx.*Ey - ny.*Ex) + Hz .* (nx.*rx + ny.*ry)) .* exp(1i*k*(rx.*rprime_x + ry.*rprime_y));
    %any NaNs here will kill the whole calculation...
    Hz_far = factor * sum(integrand(:)) * dl;
    Farfield = Hz_far;
end

return;
