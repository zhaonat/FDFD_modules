

classdef air_core_structure < handle
   %doesn't really make sense to define a single unit cell for most cases
   % also, this should be able to produce a single cell domain
   % really just a multi-grating structure... but the key field is this
   % wall_coords, which specifies the air core location
   % should be able to insert a subclass... as the wall material
   properties
       N
       xrange
       yrange
       epsilon
       Lpml
       
       wall_materials
       wall_coords     % [[y0, yf]], y0 and yf of the air core
       y_centers
       
       %% we want to be able to efficiently determine locatio specific 
       % information about this structure, particularly coord to grid for
       % source insertion or placement...
       
       %% this function really shouldn't give a damn about the materials that are 
       % used in the wall
       
   end
   
   methods
        function [obj] = air_core_structure(xrange, yrange, N, Lpml)
           obj.N = N;
           obj.xrange = xrange;
           obj.yrange=  yrange;
           obj.epsilon =ones(N); 
           %how do we deal with the pml when constructing grids?
           obj.Lpml = Lpml;
           
        end
        
        % the air core structure will work as follows
        % we give it the correct wall epsilon
        % note that wall_epsilon is created before hadn with r_start r_end
        % in mind
        % it seems weird... that me, the user has to specify the entire
        % part of the grid that contains the wall before even specifying
        % this obj?       
        
        %% if I do this...then air core thickness, location in grid is already
        %% FIXED, 
        function [] = add_wall(obj, r_start, r_end, wall_epsilon)
            %helper function
            [rnx, rny] = coord_to_grid(r_start, obj.N, obj.xrange, obj.yrange);
            [rnxf, rnyf] = coord_to_grid(r_end, obj.N, obj.xrange, obj.yrange);
            
            obj.epsilon(rnx:rnxf, rny:rnyf) = wall_epsilon;
            %store air coord coordinates
            obj.wall_materials{length(obj.wall_materials)+1} = wall_epsilon;
            obj.wall_coords{length(obj.wall_coords)+1} = {r_start, r_end};
        end
        
        
        function [] = build_core(obj, r_cell, wall_materials)
            obj.y_centers = [obj.y_centers, y_center];
            % wall_materials : 2x1 cell, each containing an epsilon subgrid
            % r_cell         : 2x1 cell, each containing an r_start, r_end.
            obj.add_wall(r_start, r_end, wall_materials{1}) 
            obj.add_wall(r_start, r_end, wall_materials{2})
            
        end
        
        function []= build_multi_core(y_centers)
            % should use the build_core function...but potential issue with 
            % shared walls
            for y_center = y_centers
               build_core
            end
            
        end
        
        function [] = build_uniform_wall_core(obj,y_center, core_thickness, wall_thickness, wall_materials)
           % we can easily define a simple function which has a uniform
           % y_center, air_core_thickness, wall_thicknesses
           obj.y_centers = [obj.y_centers, y_center];
           obj.wall_materials = wall_materials
           %wall boundaries touching the core
           y_lower_wall = y_center-core_thickness/2;
           y_upper_wall = y_center+core_thickness/2;
           
           %outer wall boundaries;
           y_lower_wall2 =y_lower_wall-wall_thickness;
           y_upper_wall2 = y_upper_wall+wall_thickness;
           
           % convention for wall coords {xcoords, ycoords} vs {coord1,
           % coord2}....ugh...
           obj.wall_coords{1} = {obj.xrange, [y_lower_wall, y_lower_wall2]};
           obj.wall_coords{2} = {obj.xrange, [y_upper_wall, y_upper_wall2]};
           
           % convert this to grid points
           %ny_center = coord_to_grid([0, y_center], obj.N, obj.xrange, obj.yrange);
           
           [~,ny_lower_wall] = coord_to_grid([0, y_lower_wall], obj.N, obj.xrange, obj.yrange);
           [~,ny_upper_wall] = coord_to_grid([0, y_upper_wall], obj.N, obj.xrange, obj.yrange);

           [~,ny_lower_wall2] = coord_to_grid([0, y_lower_wall2], obj.N, obj.xrange, obj.yrange);
           [~,ny_upper_wall2] = coord_to_grid([0, y_upper_wall2], obj.N, obj.xrange, obj.yrange);
           
           obj.epsilon(:,ny_lower_wall2:ny_lower_wall) = wall_materials(1);
           obj.epsilon(:,ny_upper_wall:ny_upper_wall2) = wall_materials(2);
           
        end
        
        
      
   end
end