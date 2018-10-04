
%% how to generate a curved geometry on a rectangular grid...
% i.e. a curved hybrid grating

%use different fractions of a circle seems the easiest solution
% so arcs

%input
% arc_angle; specified in radians since sine and cosine are in radians in
% matlab
% inner_radius
% outer radius
% circle center: we'll assume 0

function [eps] = ...
    curved_stripe(eps,N,xrange, yrange, inner_radius, outer_radius, delta_arc, eps1, eps2)
    
    %% ONLY DOES A FULL CIRCLE RIGHT NOW

    Nx = N(1); Ny = N(2);

    x = linspace(xrange(1), xrange(2), Nx);
    y = linspace(yrange(1), yrange(2), Ny);

    [X,Y] = meshgrid(x,y);
    X = X.'; Y = Y.';
    R = sqrt(X.^2+Y.^2);
    THETA = acos(X./R); %only gives values between 0 and pi...
    
    %% need a special operation since acos only goes from 0 to pi
    siphoned_half = flipud(THETA(:,ceil(Ny/2):end));
    THETA(:,ceil(Ny/2):end) =siphoned_half+pi;
    %%
    
%     test(X.^2+Y.^2 > inner_radius^2 & X.^2+Y.^2 < outer_radius^2)=12;
%     test(THETA >theta_min & THETA < theta_max) = 1;

    ring_grid_boolean=X.^2+Y.^2 > inner_radius^2 & X.^2+Y.^2 < outer_radius^2;
    %how do we enforce arc angle? need an arc grid

%     delta_arc = 18*pi/180;

    arc_strips = 0:delta_arc:2*pi;
    stripes = repmat([eps1;eps2], round(length(arc_strips)/2));

    for c = 1:length(arc_strips)-1
       min_theta = arc_strips(c);
       if(c == length(arc_strips))
          max_theta = pi/2;
       else
          max_theta = arc_strips(c+1);
       end
       eps(THETA >min_theta & THETA < max_theta & ring_grid_boolean) =stripes(c);
    end

end


% Nx = 100; Ny = 200;
% xrange = [-2,2];
% yrange = [-2,2];
% inner_radius = 1; outer_radius = 2;
% arc_angle = pi/2;
% theta_min = 0; theta_max = theta_min+arc_angle;
% x = linspace(xrange(1), xrange(2), Nx);
% y = linspace(yrange(1), yrange(2), Ny);
% %cos(theta) = x/radius; %sin(theta) = y/radius;
% %theta = arccos(x/radius);
% 
% test = ones(Nx,Ny);
% [X,Y] = meshgrid(x,y);
% X = X.'; Y = Y.';
% R = sqrt(X.^2+Y.^2);
% THETA = acos(X./R); %only gives values between 0 and pi...
% siphoned_half = flipud(THETA(:,ceil(Ny/2):end));
% THETA(:,ceil(Ny/2):end) =siphoned_half+pi;
% test(X.^2+Y.^2 > inner_radius^2 & X.^2+Y.^2 < outer_radius^2)=12;
% test(THETA >theta_min & THETA < theta_max) = 1;
% 
% figure(); imagesc(test);
% ring_grid_boolean=X.^2+Y.^2 > inner_radius^2 & X.^2+Y.^2 < outer_radius^2;
% %how do we enforce arc angle? need an arc grid
% 
% 
% %  NOW HOW DO WE CREATE A STRIPED CIRCLE?
% %  WE COULD DO IT IN WEDGES, BUT THAT MEANS INNER ARC LENGTH AND OUTER ARC
% %  LENGTH WILL DIFFER... WHICH IS FINE
% 
% delta_arc = 18*pi/180;
% 
% arc_strips = 0:delta_arc:2*pi;
% stripes = repmat([2;16], 11);
% figure();
% 
% for c = 1:length(arc_strips)-1
%    min_theta = arc_strips(c)
%    if(c == length(arc_strips))
%       max_theta = pi/2;
%    else
%       max_theta = arc_strips(c+1)
%    end
%    test(THETA >min_theta & THETA < max_theta & ring_grid_boolean) =stripes(c);
%    figure();
%    imagesc(test)
%    drawnow();
% end
