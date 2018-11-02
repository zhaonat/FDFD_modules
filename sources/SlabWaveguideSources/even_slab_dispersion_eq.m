
%% function which solves the transcendental equation for a waveguide
function [fun] = even_slab_dispersion_eq(d,omega, e_slab, e_out,mode)
    % orientation of the grid is x for propagation direction(translationally invariant)
    %z for the transverse direction
    % only solves the even modes
    
    mu0 = 4*pi*1e-7; % EVERYTHING IS SI UNITS
    
    kz_max = sqrt(omega^2*mu0*(e_slab-e_out));
    
    if(mode == 'TE')
        fun =  @(kz)tan(kz*d) -(sqrt((kz_max*d)^2-(kz*d)^2))/(kz*d);
        %tan(kzd) = alpha_z/kz = sqrt(kz_max^2-kz^2)/kz
    elseif(mode == 'TM')
        fun =  @(kz)tan(kz*d) -(e_slab/e_out)*(sqrt((kz_max*d)^2-(kz*d)^2))/(kz*d);

    else
        disp('use TE or TM')
    end
    
end