%% function generates the stretched coordinate factors for Dws
%% these apply only on boundaries

function sfactor_array = create_sfactor_mine(wrange, s, omega, eps0, mu0, ...
    Nw, Nw_pml, lnR, m)
    %% Input Parameters
    % wrange: [wmin wmax], range of domain in w-direction including PML..so
    % this is the original Nx plus 2*(pml width)
    % s: 'b' or 'f', indicating whether s-factor is for Dwb or Dwf
    % omega: angular frequency
    % eps0: vacuum permittivity
    % mu0: vacuum permeability
    % Nw: number of cells in w-direction
    % Nw_pml: number of cells in PML

    if(nargin <8)
        lnR = -16;  % R: target reflection coefficient for normal incidence
    end
    if(nargin <9)
        m = 4;% degree of polynomial grading
    end
    %% Output Parameter and Notes
    % sfactor_array: 1D array with Nw elements containing PML s-factors for Dws
    % SCALING OF THE PML S FACTOR MATTERS HUGELY IN WHETHER IT WORKS OR NOT
    %% create w coordinates
    wmin = wrange(1); wmax = wrange(2);
    w_array = linspace(wmin,wmax,Nw+1);

    right_pml = Nw_pml; %isolate points
    left_pml = length(w_array) - Nw_pml;
    %where pml will be placed
    right_pml_coord = w_array(right_pml+1);
    left_pml_coord = w_array(left_pml);
    PML_thickness = abs(wrange - [right_pml_coord, left_pml_coord]);


    wn = [right_pml,left_pml]; %coordinates where the pml starts
    eta0 = sqrt(mu0/eps0);  % vacuum impedance

    s_max = -(m+1)*lnR/(2*eta0)/PML_thickness(1);

    %% check s and create pml coordinates accordingly
    if s == 'b'
        w_array = w_array(1:end-1);
    elseif s == 'f'
        w_array = (w_array(1:end-1) + w_array(2:end))/2; %in the forward case, we are evaluating the pml
        %impedance at the midpoint
    end

    %% Finally create the s-factor array
    sfactor_array = ones(Nw,1);
    % populate the left side
    for i = 1:wn(1)
        sfactor_array(i) = 1 - 1i*s_max*((w_array(wn(1))-w_array(i))/PML_thickness(1))^m/(omega*eps0);
    end
    % populate the right side
    for i = 0:length(w_array)-wn(2)
       index = length(w_array)-i;
       sfactor_array(index) = 1 - 1i*s_max*((w_array(index)-w_array(wn(2)))/PML_thickness(1))^m/(omega*eps0);
    end


end
