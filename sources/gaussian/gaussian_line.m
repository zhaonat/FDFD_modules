
function [source] = gaussian_line(xspace, Nx, sigma)
    
    source = exp(-xspace.^2/(2*sigma)).*ones(1,Nx)/Nx;

end