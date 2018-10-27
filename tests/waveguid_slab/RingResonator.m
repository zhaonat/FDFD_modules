%% function which creates ring_resonator


function eps_ring = RingResonator(eps_r, inner_rad, outer_rad, epsilon)
    N = size(eps_r)
    xc = (round(N(1)/2)); yc = (round(N(2)/2));
    x = -xc+1:xc; y = -yc+1:yc;
    [X,Y] = meshgrid(x,y);
    eps_ring = ((X.^2+Y.^2)<outer_rad^2);
    eps_inner = ((X.^2+Y.^2)<inner_rad^2);
    eps_ring = eps_ring-eps_inner;
    eps_ring(eps_ring == 1) = epsilon;
    eps_ring(eps_ring == 0) = 1;

end