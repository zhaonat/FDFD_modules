
function [eps] = parabolic_mirror_2D(N, xrange, yrange, alpha, diel)

    %y = 4ax^2 
    % parabolas are defined by a directrix and a focus
    x = linspace(xrange(1), xrange(2), N(1));
    y = linspace(yrange(1), yrange(2), N(2));
    eps = ones(N);
    [X,Y] = meshgrid(x,y);
    eps(Y < 4*alpha*X.^2 < 1e-3) =diel;
    
    
end