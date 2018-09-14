
function [Hz_tfsf, Ex_tfsf, Ey_tfsf, source_fields_cell] =...
    tfsf_fields_TE(L0, N, wvlen,xrange, yrange, tfsf_mask, A,  Tepx, Tepy, Mz, Npml)   

    eps0 = 8.85*10^-12*L0;
    mu0 = 4*pi*10^-7*L0; 
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 

    omega = 2*pi*c0/wvlen;
    M = prod(N); Nx = N(1); Ny = N(2);
    Q = spdiags(tfsf_mask(:), 0, M,M);
    
    epsilon_uniform = ones(N);
    [A0, ~, omega, ~, ~, Dxb, Dyb] =...
        solveTE_matrices(L0, wvlen, xrange, yrange, epsilon_uniform, Mz, Npml);

    b = sparse(Mz(:));
    tic
    fsrc = A0\b; %this has no particular localization...
    toc    

    % these are the fields with no structure in it
    source_fields_Hz = reshape(fsrc,Nx,Ny);
    source_fields_Ex = reshape(1/(1i*omega) * Tepy^-1 * Dyb *fsrc, Nx,Ny);
    source_fields_Ey = reshape(1/(1i*omega) * Tepx^-1 * (-Dxb * fsrc),Nx,Ny);

    source_fields_cell = {source_fields_Hz,source_fields_Ex,source_fields_Ey};
  

    bsrc = (Q*A-A*Q)*fsrc; %Fsrc is not the same as b... which is interesting...
    tic
    hz_tfsf = A\bsrc;
    toc
    
    Hz_tfsf = reshape(hz_tfsf,Nx,Ny);
    Ex_tfsf = reshape(1/(1i*omega) * Tepy^-1 * Dyb *hz_tfsf, Nx,Ny);
    Ey_tfsf = reshape(1/(1i*omega) * Tepx^-1 * (-Dxb * hz_tfsf),Nx,Ny);

    
    
end