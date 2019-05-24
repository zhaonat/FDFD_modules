function sfactor_array = create_sfactor(wrange, s, omega, eps0, mu0, Nw, Nw_pml, lnR, m)
%% Input Parameters
    % wrange: [wmin wmax], range of domain in w-direction including PML
    % s: 'b' or 'f', indicating whether s-factor is for Dwb or Dwf
    % omega: angular frequency
    % eps0: vacuum permittivity
    % mu0: vacuum permeability
    % Nw: number of cells in w-direction
    % Nw_pml: number of cells in PML
    if(nargin <8)
        lnR = -12;  % R: target reflection coefficient for normal incidence
    end
    if(nargin <9)
        m = 3.5;% degree of polynomial grading
    end
    %% Output Parameter
    % sfactor_array: 1D array with Nw elements containing PML s-factors for Dws

    eta0 = sqrt(mu0/eps0);  % vacuum impedance
%     m = 3.5;  % degree of polynomial grading
%     lnR = -12;  % R: target reflection coefficient for normal incidence
    %prev values m = 4; lnr=-16;
    w_array = linspace(wrange(1), wrange(2), Nw+1);

    loc_pml = [w_array(1 + Nw_pml), w_array(end - Nw_pml)]; % specifies where the pml begins on each side
    d_pml = abs(wrange - loc_pml); % pml thickness

    %% what happens when we have a 0 in the denominator
    sigma_max = -(m+1)*lnR/(2*eta0) ./ d_pml; %usually the pml is the same thickness on both sides

    %% forward or backward...idk what this is requiring
    if s == 'b'
        ws = w_array(1:end-1);
    else  % s == 'f'
        assert(s == 'f');
        ws = (w_array(1:end-1) + w_array(2:end)) / 2;
    end

    ind_pml = {ws < loc_pml(1), ws > loc_pml(2)};  % {} means cell array

    sfactor_array = ones(1, Nw);
    for n = 1:2
        sfactor = @(L) 1 - 1i * sigma_max(n)/(omega*eps0) * (L/d_pml(n)).^m;
        sfactor_array(ind_pml{n}) = sfactor(abs(loc_pml(n) - ws(ind_pml{n})));
    end

end