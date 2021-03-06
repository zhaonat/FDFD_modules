
%% this class has the ability to insert multiple gratings into itself
%% which is sort of inconsistent from the class point of view

classdef periodic_grating < handle
   %doesn't really make sense to define a single unit cell for most cases
   % also, this should be able to produce a single cell domain
   % structure also assumes that the grating periodicity parallel x direction
   properties
       N
       xrange
       yrange
       epsilon
       Lpml
       
       %% Do we want the ability to hold multiple grating specs here or in the air core
       %% it has to be here since air_core inherits the periodic grating
       %% however, it's weird that periodic grating can add multiple walls by itself
       grating_properties = cell(0);
       
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
        
        % all the subsequent functions assume a horizontally oriented
        % grating array, which is somewhat restricted...
        % this is the main function to use to generate a periodic
        % grating...
        function [] = add_grating_array(obj, num_cells, lattice_constant, ...
                thickness, epsilon_array, fill_factor, y_center)
            
            assert(length(epsilon_array) == 2, 'need two distinct eps_r') 
            assert(num_cells*lattice_constant <= diff(obj.xrange), 'num cells is too large or lattice constant too large');

            wall_properties = struct;
            wall_properties.eps1 = epsilon_array(1);
            wall_properties.eps2 = epsilon_array(2);
            wall_properties.lattice_constant = lattice_constant;
            wall_properties.num_cells = num_cells;
            wall_properties.y_pos = y_center;
            wall_properties.wall_coords = [y_center-thickness/2, y_center+thickness/2];
            wall_properties.grating_thickness = thickness;
            
                                    
            for i = 0:num_cells-1
               x_center = obj.xrange(1)+(i)*lattice_constant+(lattice_constant/2);
               wall_properties=obj.add_horizontal_grating(epsilon_array(1), epsilon_array(2), ...
                    fill_factor, thickness, y_center, x_center, lattice_constant, wall_properties);
            end
            obj.grating_properties{length(obj.grating_properties)+1} = wall_properties;
            obj.remove_redundant_properties();
            
            %sizing constraints
            if(~isequal(size(obj.epsilon), obj.N))
                obj.epsilon = obj.epsilon(1:obj.N(1), 1:obj.N(2));
            end
            
        end
        
        %% helper function adds a single unit cell
        % we can hijack this function by externally transposing all the
        % inputs...
        function [wall_properties] = ...
                add_horizontal_grating(obj, epsilon_diel, epsilon_metal, ...
            fill_factor, thickness, y_center, x_center, lattice_constant, wall_properties)
        
            L = [diff(obj.xrange), diff(obj.yrange)];
            lattice_n = (lattice_constant/L(1))*(obj.N(1));

            [nxcenter, nycenter] =...
                coord_to_grid([x_center, y_center],obj.N, obj.xrange,obj.yrange)
            metal_size = fill_factor*lattice_n;
            halfx = floor(metal_size/2); %floor is the correct option
            lat_half = floor(lattice_n/2);

            grating_size = round(obj.N(2)*(thickness/L(2)));
            halfy = round(grating_size/2);
            
            %check if lat_half = nxcenter
            if(lat_half == nxcenter)
               nxcenter = nxcenter+1; 
            end
            wall_properties.nygrid = [nycenter-halfy,nycenter+halfy];
            %% eps2 part
            obj.epsilon(nxcenter-halfx: nxcenter+halfx, nycenter-halfy:nycenter+halfy) = epsilon_metal;
            %% eps1 parts
            obj.epsilon(nxcenter-lat_half:nxcenter-halfx, nycenter-halfy:nycenter+halfy) = epsilon_diel;
            obj.epsilon(nxcenter+halfx:nxcenter+lat_half, nycenter-halfy:nycenter+halfy) = epsilon_diel;
% 
%             nybounds = [nycenter-halfy, nycenter+halfy];
%             metal_bounds_nx = [nxcenter-halfx, nxcenter+halfx];

        end   
        
        function remove_redundant_properties(obj)
            for i = 1:length(obj.grating_properties)-1
                if(isequal(obj.grating_properties{i}, obj.grating_properties{i+1}))
                   obj.grating_properties{i} = [];
                end
                obj.grating_properties=obj.grating_properties(~cellfun('isempty',obj.grating_properties));
            end

            
        end
        
   end
end