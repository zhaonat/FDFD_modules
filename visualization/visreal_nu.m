function h = visreal_nu(array2d, xrange_array, yrange_array, grid_line_bool)

    %% Get the maximum magnitude before taking the real part.
    cmax = max(abs(array2d(:)));

    %% Attach a row and column at the ends.  (Not visualized, but required by pcolor().)
    array2d = real(array2d);

    %% Create the matrices for x-locations (X), y-locations (Y), and color (C).
    [X, Y] = meshgrid(xrange_array, yrange_array);
    C = permute(array2d, [2 1]);

    %% Draw with pcolor().
    h = pcolor(X, Y, C);

    %% Make the figure look better.
    if(grid_line_bool==0)
        set(h, 'EdgeColor', 'none');
    end
    set(gca, 'TickDir', 'out');
    axis image;

    caxis([-cmax, cmax]);
    colormap('b2r')
    colorbar;

end
