
function eps_proto = createCube(eps_vac, edgeLength, epdiel)
    
    dim = size(eps_vac);
    Nx = dim(1); Ny = dim(2); Nz = dim(3);
    midptx = ceil(Nx/2); midpty = ceil(Ny/2); midptz = ceil(Nz/2);
    eps_vac(midptx-ceil(edgeLength/2)+1:midptx+ceil(edgeLength/2), ...
        midpty-ceil(edgeLength/2)+1:midpty+ceil(edgeLength/2), ...
        midptz-ceil(edgeLength/2)+1:midptz+ceil(edgeLength/2)) = epdiel;
    eps_proto = eps_vac;

end