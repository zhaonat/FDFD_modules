

classdef periodic_grating < handle
   %doesn't really make sense to define a single unit cell for most cases
   % also, this should be able to produce a single cell domain
    
   properties
       N
       xrange
       yrange
       epsilon
       Lpml
       
       num_cells
       lattice_constant
       eps1 % epsilon corresponding to 1-fill_fraction
       eps2 % epsilon corresponding to fill_fraction
       thickness
       x_pos
       y_pos       % center line of the periodic array
       fill_factor % defines the ratio of eps2 to eps1
       
       %% we want to be able to efficiently determine locatio specific 
       % information about this structure, particularly coord to grid for
       % source insertion or placement...
       
   end
   
   methods
        function [obj] = periodic_grating(xrange, yrange, N, Lpml)
           obj.N = N;
           obj.xrange = xrange;
           obj.yrange=  yrange;
           obj.epsilon =ones(N); 
           
           %how do we deal with the pml when constructing grids?
           obj.Lpml = Lpml;
           
        end
        
        % this is the main function to use to generate a periodic
        % grating...
        function [] = add_grating_array(obj, num_cells, lattice_constant, ...
                thickness, epsilon, fill_factor, y_center)
            obj.eps1 = epsilon(1);
            obj.eps2 = epsilon(2);
            obj.lattice_constant = lattice_constant;
            obj.num_cells = num_cells;
            obj.y_pos = y_center;
            assert(num_cells*lattice_constant <= diff(obj.xrange), 'num cells is too large or lattice constant too large');
            for i = 0:obj.num_cells-1
               x_center = obj.xrange(1)+(i)*lattice_constant+(lattice_constant/2)
               obj.add_grating(epsilon(1), epsilon(2), ...
                    fill_factor, thickness, y_center, x_center, lattice_constant);
            end
            
        end
        
        %% helper function adds a single unit cell
        function [] = ...
                add_grating(obj, epsilon_diel, epsilon_metal, ...
            fill_factor, thickness, y_center, x_center, lattice_constant)
        
            L = [diff(obj.xrange), diff(obj.yrange)];
            lattice_n = (lattice_constant/L(1))*(obj.N(1))

            [nxcenter, nycenter] =...
                coord_to_grid([x_center, y_center],obj.N, obj.xrange,obj.yrange)
            
            metal_size = fill_factor*lattice_n;
            halfx = floor(metal_size/2); %floor is the correct option
            lat_half = floor(lattice_n/2)

            grating_size = round(obj.N(2)*(thickness/L(2)));
            halfy = round(grating_size/2);

            %% eps2 part
            obj.epsilon(nxcenter-halfx: nxcenter+halfx, nycenter-halfy:nycenter+halfy) = epsilon_metal;
            %% eps1 parts
            obj.epsilon(nxcenter-lat_half:nxcenter-halfx+1, nycenter-halfy:nycenter+halfy) = epsilon_diel;
            obj.epsilon(nxcenter+halfx:nxcenter+lat_half, nycenter-halfy:nycenter+halfy) = epsilon_diel;
% 
%             nybounds = [nycenter-halfy, nycenter+halfy];
%             metal_bounds_nx = [nxcenter-halfx, nxcenter+halfx];

        end       
   end
end