

function [filtered_modes, filtered_eigs, mask] = ...
    mode_filtering(eigenmodes, eigenvals, ...
    structure_xbounds, structure_ybounds, L, Npml, pml_threshold)

    %% assumes that xrange and yrange is [-, +] and is centered on 0
    % xlim: [x0, xf] x bounds of the STRUCTURE in PHYSICAL UNITS
    % ylim: [y0, yf] y bounds of the STRUCTURE in PHYSICAL UNITS (microns
    % or whatever)
    % eigenmodes should be a cell where each cell index is a field pattern
    % or mode pattern
    if(nargin <8)
       pml_threshold = 1e-4; 
    end
    
    % mask KEYS
    % 2 = pml
    % 1 = structure
    % 0 = air;
    N = size(eigenmodes{1});
    Nx = N(1); Ny = N(2);
    Nxc = round(Nx/2); Nyc = round(Ny/2);
    x0 = structure_xbounds(1); xf = structure_xbounds(2); y0=structure_ybounds(1); yf = structure_ybounds(2);
    %convert the physical bounds to grid bounds

    Nx0 = Nxc+round((x0/L(1))*N(1))+1; Nxf = Nxc+floor((xf/L(1))*N(1));
    Ny0 = Nyc+round((y0/L(2))*N(2))+1; Nyf = Nyc+floor((yf/L(2))*N(2));

    %% get PML bounds
    x = 1:Nx; y = 1:Ny; % x and y are node grids
    [X,Y] = meshgrid(x,y);
    X = X.'; Y = Y.';
    mask = zeros(N);
    mask((X<Npml(1) | X > Nx-Npml(1)) | ...
            (Y<Npml(2) | Y> Ny - Npml(2))) = 2;
    mask(Nx0:Nxf, Ny0:Nyf) = 1;
    size(mask)
    n = length(eigenmodes);
    filtered_eigs = [];
    filtered_modes = [];
    c = 1;
    
    %% should we do an epsilon map of pml, air, and structure fields?
    
    %% execute the filters
    for i = 1:n
        
        % execute a filter which attempts to cut out modes that oscillate
        % too quickly in space (usually a problem for plasmonics)
        
        
        structure_fields = eigenmodes{i}(mask == 1);
        %get fields outside of structure
        
        air_fields = eigenmodes{i}(mask == 0);
        
        %get fields inside structure
        PML_fields = eigenmodes{i}(mask == 2);
       
        if(mean(mean(abs(PML_fields)))>pml_threshold)
            disp('pml fields too large')
            continue;   
        end

        if(mean(abs(structure_fields))> mean(abs(air_fields)))
            filtered_eigs(c) =  eigenvals(i);
            filtered_modes{c} = eigenmodes{i};
            c = c+1;
        else
            disp('too much field outside')
        end 

    end

end