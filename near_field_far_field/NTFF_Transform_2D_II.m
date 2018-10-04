function Farfield = NTFF_Transform_2D_II(E_near, H_near, norm_surf, r_near, ...
    dL, r_hat, wavelength, polarization)
    % this function should be able to handle different dx and dy as dL = [dx,
    % dy]
    % function Farfield = NTFF_Transform_2D(E, H, norm_surf, r_near, dl, r_hat, wavelength, polarization)
    % If s-polared, E = Ez [1*n], 
    %H = [Hx; Hy] [2*n], if p-polarized, E = {Ex,
    % Ey}, H = Hz; norm_surf [2*n], r_near [rx; ry], r_hat is observe direction 
    % [2*1]; 
    %dl is a scalar... but what should it be?
    %dl should actually not be a scalar...
    %H is the normalized H (Z0H). Ignore phasor exp(-jkr)/sqrt(r).

    k = 2*pi/wavelength;
    factor = exp(1i*pi/4)*sqrt(k/8/pi);
    dx = dL(1); dy = dL(2);
    dL_vec = [dx, dy, dx, dy]; %up, left, down, right
    %dL_vec = [dy, dx, dy, dx];
    % comes as a vector still
    rx = r_hat(1, :);
    ry = r_hat(2, :);

    if polarization == 's'
        integral = 0;

        for i = 1:4
            nx = norm_surf{i}(1,:);
            ny = norm_surf{i}(2,:);
            rprime_x = r_near{i}(1,:);
            rprime_y = r_near{i}(2,:);
            dl = dL_vec(i);
            Ez = E_near{i}(3,:);
            Hx = H_near{i}(1,:);
            Hy = H_near{i}(2,:);
            integrand = dl*((- (nx.*Hy - ny.*Hx) + Ez .* (nx.*rx + ny.*ry)) .* exp(1i*k*(rx.*rprime_x + ry.*rprime_y)));
            integral = integral+sum(integrand(:));
        end
        Ez_far = factor * integral;
        Farfield = Ez_far;
    elseif polarization == 'p'
        integral = 0;
        for i = 1:4
            nx = norm_surf{i}(1,:);
            ny = norm_surf{i}(2,:);
            rprime_x = r_near{i}(1,:);
            rprime_y = r_near{i}(2,:);
            dl = dL_vec(i);
            Ex = E_near{i}(1,:);
            Ey = E_near{i}(2,:);
            Hz = H_near{i}(3,:);
            integrand = dl*(((nx.*Ey - ny.*Ex) + Hz .* (nx.*rx + ny.*ry)) .* exp(1i*k*(rx.*rprime_x + ry.*rprime_y)));
            integral = integral+sum(integrand(:));
        end
        %any NaNs here will kill the whole calculation...
        Hz_far = factor * integral;
        Farfield = Hz_far;
    end

return;
