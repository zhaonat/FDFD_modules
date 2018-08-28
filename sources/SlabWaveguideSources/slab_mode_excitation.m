% function directly assumes you've already solved for the correct k and
% kappa 
function Jz_line = slab_mode_excitation(xrange, h, k, kappa, right, left)
    num_points = diff(xrange)/h;
    Jz_line = zeros(num_points,1);
    X = linspace(xrange(1), xrange(2), num_points);
    % convert r and left to num_points
    % r = xrange(1)+(r/length)*num_points;
    % left = xrange(1)+(left/length)*num_points;
    d = right-left;
    %how do we enforce continuity?
    for i = 1:num_points
        x = X(i);
        if(x > right)
            Jz_line(i) = exp(-kappa*x);
        elseif(x < left)
           Jz_line(i) = exp(kappa*x);
        else
            Jz_line(i) = cos(k*x/(2*d));
        end
    end
end

%% TEST CODE
% h = 0.1;
% xrange = [-10,10];
% num_points = diff(xrange)/h;
% Jz_line = zeros(num_points,1);
% length = diff(xrange);
% X = linspace(xrange(1), xrange(2), num_points);
% r = 2; left = -2;
% % convert r and left to num_points
% % r = xrange(1)+(r/length)*num_points;
% % left = xrange(1)+(left/length)*num_points;
% slab_width = 2;
% 
% x0 = [1,2]; %% initial guess
% omega = 1;
% [k_soln,fval] = fsolve(@(x) dielectricSlabEq(x, omega), x0)
% 
% kappa = k_soln(1);
% k = k_soln(2);
% %kappa and k satisfy a wave equation
% 
% %how do we enforce continuity?
% % solve this equation kappa = kx tan(kxd) 
% %and                  kappa^2 + k^2 = omega^2*(e1-e2)
% for i = 1:num_points
%     x = X(i);
%     if(x > r)
%         Jz_line(i) = exp(-kappa*x)
%     elseif(x < left)
%        Jz_line(i) = exp(kappa*x)
%     else
%         Jz_line(i) = cos(k*x)
%     end
% end
% 
% figure; 
% plot(X,Jz_line)