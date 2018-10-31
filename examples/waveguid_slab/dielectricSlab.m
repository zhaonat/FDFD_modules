function eps_slab = dielectricSlab(eps_r, dL, width, epsilon)
    % width is in nanometers
    N = size(eps_r);
    dx = dL(1); dy = dL(2);
    node_slab = width/dy;
    %xc = round(N(1)/2);
    yc = round(N(2)/2);
    eps_r(:,yc-round(node_slab/2):yc+round(node_slab/2),:)=epsilon;
    eps_slab = eps_r;
end