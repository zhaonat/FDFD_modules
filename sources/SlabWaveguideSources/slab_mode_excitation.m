% function directly assumes you've already solved for the correct k and
% kappa 
function Jz_line = slab_mode_excitation(xrange, h, k, kappa, right, left)
    num_points = diff(xrange)/h;
    Jz_line = zeros(num_points,1);
    X = linspace(xrange(1), xrange(2), num_points);
    % convert r and left to num_points
    % r = xrange(1)+(r/length)*num_points;
    % left = xrange(1)+(left/length)*num_points;
    %how do we enforce continuity; we have to solve for the front
    %coefficient.
    % we assume B = 1
    % field ~Acos(k_z z) or Be^(-alpha z) assume B = a;
    d = abs(right-left)/2;
    A = exp(-kappa*d)/cos(k*d);
    for i = 1:num_points
        x = X(i);
        if(x > right)
            Jz_line(i) = exp(-kappa*x);
        elseif(x < left)
           Jz_line(i) = exp(kappa*x);
        else
            Jz_line(i) = A*cos(k*x);
        end
    end
end