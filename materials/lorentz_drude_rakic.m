
function [epsilon] = lorentz_drude_rakic(lambda_scan, f0, oscillator_f, omega_p, gamma_p, omega_l_arr, gamma_l_arr)
    
    %lambda scan must be in SI UNITS
    c0 = 3e8;
    eV = 241.79893*1e12*(2*pi);
    rad_Hz = 1/(2*pi);
    Hz_eV =  4.13566553853599E-15; %1Hz = x eV
    rad_eV = rad_Hz*Hz_eV;
    
    omega_scan = 2*pi*c0./lambda_scan*rad_eV;
    big_omega_p = sqrt(f0)*omega_p;

    drude_part = 1-(big_omega_p^2)./(omega_scan.*(omega_scan+1i*gamma_p));
    
    num_resonances = length(omega_l_arr);
    lorentz_part = 0;
    for j = 1:num_resonances
        omega_j = omega_l_arr(j);
        gamma_j = gamma_l_arr(j);
        lorentz_part = lorentz_part + ...
            oscillator_f(j)*omega_p^2./(omega_j^2 - omega_scan.^2 -1i*gamma_j*omega_scan);
    end

   epsilon = drude_part + lorentz_part;

end