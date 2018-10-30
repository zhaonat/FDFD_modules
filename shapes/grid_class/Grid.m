
classdef Grid
    properties
       N
       xrange
       yrange
       Npml
       dL
    end
    methods
      function obj = Grid(xrange, yrange, N)
        obj.xrange = xrange;
        obj.yrange = yrange;
        obj.N = N;
        obj.dL(1) = diff(xrange)/N(1);
        obj.dL(2) = diff(yrange)/N(2);
      end
      
      % functions to convert coordinates to ncoords and ncoords to (x,y)
      function [nx, ny] = coord_to_node(obj, x, y)
        Lx = diff(obj.xrange); Ly = diff(obj.yrange);
        xfrac = (x-obj.xrange(1))/Lx;
        yfrac = (y-obj.yrange(1))/Ly;
        nx = round(xfrac*obj.N(1));
        ny = round(yfrac*obj.N(2));
      end
      
      function [x,y] = node_to_coord(obj, nx, ny)
          xnfrac = nx/obj.N(1);
          ynfrac = ny/obj.N(2);
          x = obj.xrange(1)+xnfrac*diff(obj.xrange);
          y = obj.yrange(1)+ynfrac*diff(obj.yrange);
          
      end
      
      function [xrange, yrange, N, dL, Lpml] = domain_with_pml(obj)
        %% Input parameters
        % xrange: [xmin xmax], range of domain in x-direction without PML
        % yrange: [ymin ymax], range of domain in y-direction without PML
        % N: [Nx Ny], numbers of cells in x- and y-directions without PML
        % Npml: [Nx_pml Ny_pml], numbers of cells in the x- and y-normal PML

        L = [diff(obj.xrange) diff(obj.yrange)];  % [Lx Ly]
        dL = L./obj.N;  % [dx dy]

        Lpml = obj.Npml .* dL;  % [Lx_pml, Ly_pml]
        xrange = obj.xrange + [-1 1] * Lpml(1);  % [xmin xmax] is updated
        yrange = obj.yrange + [-1 1] * Lpml(2);  % [ymin ymax] is updated

        N = obj.N + 2*obj.Npml;  % [Nx Ny] is updated
      end
      
      
      
    end
      
      
end
