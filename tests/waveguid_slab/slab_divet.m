function eps_divet = slab_divet(eps_vac,dL, epsilon, ...
    divet_width, divet_height, slab_width)
    N = size(eps_vac);
    xc= round(N(1)/2); yc = round(N(2)/2);
    eps_slab = dielectricSlab(eps_vac,dL, slab_width, epsilon );
    
    N_slab_width = slab_width/dL(1);
    N_divet_height = divet_height/dL(1);
    N_divet_width = divet_width/dL(1);

    %% divet
    eps_slab(xc-N_divet_width:xc+N_divet_width,...
        yc+round(N_slab_width/2):yc+round(N_slab_width/2)+N_divet_height) = epsilon;
    
    eps_divet = eps_slab;
end