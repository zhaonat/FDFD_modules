%% function generates the stretched coordinate factors for Dws
%% these apply only on boundaries

function sfactor_array = create_sfactor_mine(wrange, s, omega, eps0, mu0, Nw, Nw_pml)
    %% Input Parameters
    % wrange: [wmin wmax], range of domain in w-direction including PML..so
    % this is the original Nx plus 2*(pml width)
    % s: 'b' or 'f', indicating whether s-factor is for Dwb or Dwf
    % omega: angular frequency
    % eps0: vacuum permittivity
    % mu0: vacuum permeability
    % Nw: number of cells in w-direction
    % Nw_pml: number of cells in PML

    %% Output Parameter and Notes
    % sfactor_array: 1D array with Nw elements containing PML s-factors for Dws
    % SCALING OF THE PML S FACTOR MATTERS HUGELY IN WHETHER IT WORKS OR NOT
    %% create w coordinates
    wmin = wrange(1); wmax = wrange(2);
    wcoord = linspace(wmin,wmax,Nw+1);

    right_pml = Nw_pml; %isolate points
    left_pml = length(wcoord) - Nw_pml;
    %where pml will be placed
    right_pml_coord = wcoord(right_pml+1);
    left_pml_coord = wcoord(left_pml);
    PML_thickness = abs(wrange - [right_pml_coord left_pml_coord]);


    wn = [right_pml,left_pml]; %coordinates where the pml starts
    eta0 = sqrt(mu0/eps0);  % vacuum impedance
    m = 4;  % degree of polynomial grading
    lnR = -16;  % R: target reflection coefficient for normal incidence

    s_max = -(m+1)*lnR/(2*eta0)/PML_thickness(1);

    %% check s and create pml coordinates accordingly
    if s == 'b'
        wcoord = wcoord(1:end-1);
    elseif s == 'f'
        wcoord = (wcoord(1:end-1) + wcoord(2:end))/2; %in the forward case, we are evaluating the pml
        %impedance at the midpoint
    end

    %% Finally create the s-factor array
    sfactor_array = ones(Nw,1);
    for i = 1:wn(1)
        sfactor_array(i) = 1 - 1i*s_max*((wcoord(wn(1))-wcoord(i))/PML_thickness(1))^m/(omega*eps0);
    end
    for i = 0:length(wcoord)-wn(2)
       index = length(wcoord)-i;
       sfactor_array(index) = 1 - 1i*s_max*((wcoord(index)-wcoord(wn(2)))/PML_thickness(1))^m/(omega*eps0);
    end

end
