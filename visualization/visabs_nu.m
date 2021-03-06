function visabs_nm(array2d, xs, ys)

%% Attach a row and column at the ends.  
% The extra row and column are not visualized, but required by pcolor().
% array2d = log(abs(array2d)); 
size(array2d)

array2d = abs(array2d);
[Nx, Ny] = size(array2d);
% array2d = [array2d, array2d(:,1)];
% array2d = [array2d; array2d(1,:)];

%% Create the matrices for x-locations (X), y-locations (Y), and color (C).
% xs = linspace(xrange(1), xrange(2), Nx+1);
% ys = linspace(yrange(1), yrange(2), Ny+1);
% [X, Y] = meshgrid(xs, ys);
% xs = cumsum(dx_vec); 
% ys = cumsum(dy_vec); 
[X, Y] = meshgrid(xs, ys); 


size(array2d)

C = permute(array2d, [2 1]);

%% Draw with pcolor().
h = pcolor(X, Y, C);

%% Make the figure look better.
set(h, 'EdgeColor', 'none');
set(gca, 'TickDir', 'out');
axis image;

colormap('hot')
colorbar;
