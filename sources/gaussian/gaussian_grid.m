
function [source] = gaussian_grid(N,Npml, xrange, yrange, sigma,Kx)
    %N0 = N-Npml;
    x = linspace(xrange(1), xrange(2), N(1));
    y = linspace(yrange(1), yrange(2), N(2));
    [X,Y] = meshgrid(x,y);
    X = X.'; Y = Y.';
    source = exp(-(Y.^2)/sigma^2).*exp(-1i*Kx*X);
    source(1:Npml(1),:) = 0;
    source(N(1)-Npml(1):end,:) = 0;   
    source(:,1:Npml(2)) = 0;
    source(:,N(2)-Npml(2):end) = 0;
end