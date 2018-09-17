
function [eps] = parabolic_mirror_2D(eps_multi,  xrange_array, yrange_array, alpha, diel)

    %y = 4ax^2 
    % parabolas are defined by a directrix and a focus
    % for use with our hybrid grating, we would actually like to construct
    % a bulk rectangular dielectric structure and then chop off the part
    % which is irrelevant
    
%     x = linspace(xrange(1), xrange(2), N(1));
%     y = linspace(yrange(1), yrange(2), N(2));
	%N = size(eps_multi);
    [X,Y] = meshgrid(xrange_array,yrange_array);
    eps(Y < 4*alpha*X.^2 < 1e-3) =diel;
    
    
end