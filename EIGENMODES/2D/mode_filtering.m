

function [filtered_modes, filtered_eigs] = ...
    mode_filtering(eigenmodes, omega_eigs,  eps, xlim, ylim, L, Npml)

    %
    % xlim: [x0, xf] xbounds 
    % ylim: [y0, yf] y bounds
    
    % mask code
    % 2 = pml
    %1 = structure
    % 0 = air;
    N = size(eps);
    Nx = N(1); Ny = N(2);
    xc = round(Nx/2); yc = round(Ny/2);
    x0 = xlim(1); xf = xlim(2); y0=ylim(1); yf = ylim(2);
    %convert the physical bounds to grid bounds

    Nx0 = xc+round((x0/L(1))*N(1))+1; Nxf = xc+round((xf/L(1))*N(1));
    Ny0 = yc+round((y0/L(2))*N(2))+1; Nyf = yc+round((yf/L(2))*N(2));
    
    %% get PML bounds
    x = 1:Nx; y = 1:Ny;
    [X,Y] = meshgrid(x,y);
    
    mask = zeros(N);
    mask((X<Npml(1) | X > Nx-Npml(1)) | ...
            (Y<Npml(2) | Y> Ny - Npml(2))) = 2;
    mask(Nx0:Nxf, Ny0:Nyf) = 1;
    n = length(eigenmodes);
    filtered_eigs = [];
    filtered_modes = [];
    c = 1;
    
    %% should we do an epsilon map of pml, air, and structure fields?
    
    %% execute the filters
    for i = 1:n
        structure_fields = eigenmodes{i}(mask == 1);
        %get fields outside of structure
        
        air_fields = eigenmodes{i}(mask == 0);
        
        %get fields inside structure
        PML_fields = eigenmodes{i}(mask == 2);
       
        if(mean(mean(abs(PML_fields)))>1e-1)
            continue;   
        end

        if(mean(abs(structure_fields))> mean(abs(air_fields)))
            filtered_eigs = [filtered_eigs, omega_eigs(i)];
            filtered_modes{c} = eigenmodes{i};
            c = c+1;
        end 

    end

end