function [r_near, norm_surf, E_near, H_near] =...
    Get_Nearfield_2D(box_bound, E, H, xrange, yrange, dL, side_names, polarization)
% function [r_near, norm_surf, E_near, H_near] = Get_Nearfield_2D(box_bound, E, H, xrange, yrange, grid_sample, side_names, polarization)
%% input
% by zhexin zhao
% it appears that the fields in the PML cannot be interpolated?
% box_bound = [xmin, xmax; ymin, ymax]; %in the physical units
% E and H are results from FDFD and are in their flattened form
% xrange, yrange, are domain size in simulation [min, max]
% side_names = {'up', 'right', 'left', 'down'} with the sides included in
% NTFF 
% polarization = 'p' or 's'. If 'p', H is a scalar, E is a cell of 2 elems.
% If 's', E is a scalar, H is a cell of 2 elements

%% output
% r_near: an array of near field positions, each colume is [rx; ry]
% norm_surf: an array of surface norm vectors, each colume is [nx, ny]
% E_near, H_near: near fields, [Ex; Ey; Ez]

%% get data size N = [Ny, Nx];
if polarization == 's'
    N = size(E);
elseif polarization == 'p'
    N = size(H);
end

%% obtain side flags
% flag == 0: not include, flag == 1: included
% side_names are places where we want to get the near-field...so up...
up_flag = 0;
left_flag = 0;
right_flag = 0;
down_flag = 0;
for ind = 1:size(side_names, 2)
    if strcmp(side_names{ind}, 'up')
        up_flag = 1;
    elseif strcmp(side_names{ind}, 'left')
        left_flag = 1;
    elseif strcmp(side_names{ind}, 'right')
        right_flag = 1;
    elseif strcmp( side_names{ind}, 'down')
        down_flag = 1;
    else
        disp('Error: Wrong input in side_names!')
    end
end

%% get dx, dy
dX = (xrange(2) - xrange(1)) / N(2);
dY = (yrange(2) - yrange(1)) / N(1);
sx = linspace(xrange(1), xrange(2)-dX, N(2));
sy = linspace(yrange(1), yrange(2)-dY, N(1));

[sx_prime, sy_prime] = meshgrid(sx, sy);
[sx_dual, sy_dual] = meshgrid(sx+dX/2, sy+dY/2); %averaging needed
dx = dL(1);
dy = dL(2);

%% get nearfield sampling positions
xmin = box_bound(1,1);
xmax = box_bound(1,2);
ymin = box_bound(2,1);
ymax = box_bound(2,2);

% up
x_up = (xmin+dx/2):dx:(xmax - dx/2);
y_up = ymax*ones(size(x_up));
% left
y_left = (ymin+ dy/2):dy:(ymax - dy/2);
x_left = xmin*ones(size(y_left));
% right
y_right = (ymin + dy/2):dy:(ymax - dy/2);
x_right = xmax*ones(size(y_right));
% down
x_down = (xmin + dx/2):dx:(xmax - dx/2);
y_down = ymin*ones(size(x_down));

% create a mask for different sizes
%[inverse(up), inverse(left), down, right]
mask_flag = [up_flag * ones(size(x_up)), left_flag * ones(size(y_left)), ...
            down_flag * ones(size(x_down)), right_flag * ones(size(y_right))];
%mask_flag is 1x4*
%% r_near contains the coordinates of the bounding box used to get the near-fields
r_near = [[x_up(end:(-1):1); y_up], ...
    [x_left; y_left(end:(-1):1)], ...
    [x_down; y_down],...
    [x_right; y_right]];

% we want to separate r_near into four separate containers

% r_near{1} = [x_up(end:(-1):1); y_up];     %up
% r_near{2} = [x_left; y_left(end:(-1):1)]; %left
% r_near{3} = [x_down; y_down];             %down
% r_near{4} = [x_right; y_right];           %right

norm_surf = [[zeros(size(x_up)); ones(size(x_up))], ...  % up: (0;1)
    [-1*ones(size(y_left)); zeros(size(y_left))],...    % left: (-1;0)
    [zeros(size(x_down)); -1*ones(size(x_down))], ...   % down: (0;-1)
    [ones(size(y_right)); zeros(size(y_right))]];       % right: (1;0)

%% get field values at sampling positions
if polarization == 's'
    Ez = interp2(sx_prime, sy_prime, E, r_near(1,:), r_near(2,:));      % Ez field at sampling points
    Hx = interp2(sx_prime, sy_dual, H{1}, r_near(1,:), r_near(2,:));   % Hx
    Hy = interp2(sx_dual, sy_prime, H{2}, r_near(1,:), r_near(2,:));   % Hy
    
    %% WHY do we construct three rows for the near field vector?
    % it's because it's for one for each component...apparently.
    E_near = [zeros(size(Ez)); zeros(size(Ez)); mask_flag.* Ez];
    % see Hx, Hy, Hz...which is 0 for 's'
    H_near = [mask_flag.* Hx; mask_flag.* Hy; zeros(size(Hx))];
elseif polarization == 'p'
    % r_near(1,:) = Xq; r_near(2,:) = Yq;
    Ex = interp2(sx_dual, sy_prime, E{1}, r_near(1,:), r_near(2,:));   % Ex
    Ey = interp2(sx_prime, sy_dual, E{2}, r_near(1,:), r_near(2,:));   % Ey
    Hz = interp2(sx_dual, sy_dual, H, r_near(1,:), r_near(2,:));      % Hz
    % Ex, Ey, Ez (which is 0 for P polarization)
    E_near = [mask_flag.* Ex; mask_flag.* Ey; zeros(size(Ex))];
    H_near = [zeros(size(Hz)); zeros(size(Hz)); mask_flag.* Hz];
else
    disp('Polarization error!')
end

return;
