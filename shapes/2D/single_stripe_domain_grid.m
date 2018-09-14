
function eps = single_stripe_domain_grid(N,L,Delta_stripe, eps_1, eps_2)
    %N: grid points [Nx,Ny]
    %L: [Lx,Ly]; physical dimensions of grid
    % Delta_Stripe [thick_x, thick_y]; size of the stripe
    % eps_1; background dielectric
    % eps_2: stripe dielectric
    dL = L./N;
    
    eps = ones(N);
    ox = floor(N(1)/2);
    oy = floor(N(2)/2); % origin
    
    Delta_x = Delta_stripe(1); Delta_y = Delta_stripe(2);
    halfx = floor(Delta_x/dL(1))/2;
    halfy = floor(Delta_y/dL(2))/2;
    
    eps(ox-halfx+1:ox+halfx,:) = eps_1;
    eps(ox-halfx+1:ox+halfx, oy-halfy+1: oy+halfy) = eps_2;


end