
function Jz_plane = PlaneWaveSource(k, dL, location_mask,N)
    x = 1:N(1); y = 1:N(2);
    kx = k(1); ky = k(2);
    hx = dL(1); hy = dL(2);
    [X,Y] = meshgrid(x,y);
    
    X(location_mask==0) = 0;
    Y(location_mask==0) = 0;

    Jz_plane = exp(1i*kx*X*hx + 1i*ky*Y*hy); 
    Jz_plane(location_mask==0) = 0;
    
end
