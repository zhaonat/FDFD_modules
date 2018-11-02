

classdef air_core_structure < handle & periodic_grating
   % this should really inherit from periodic grating
   % since this is, in its most complex form, just a set of periodic
   % gratings. However, we should not inherit it periodic grating because
   % then the field
   
   % really just a multi-grating structure... but the key field is this
   % wall_coords, which specifies the air core location
   % should be able to insert a subclass... as the wall material
   
   % all methods and properties are now inherited from periodic grating
   properties       
       core_centers_y
       wall_bounds
       
   end
   
   methods
        function [obj] = air_core_structure(xrange, yrange, N, Lpml)
             obj@periodic_grating(xrange, yrange, N,Lpml)
        end
    
        
        function []= build_multi_core(obj, y_centers,  wall_materials_cell)
            % should use the build_core function...but potential issue with 
            % shared walls
            % y_centers....these are the central coordinates of walls
            assert(length(y_centers) == length(wall_materials_cell), ...
                'number of materials must be same as number of walls')
            
            for i = 1:length(y_centers) 
                
               y_wall_center = y_centers(i);
               wall_properties = wall_materials_cell{i};
               [num_cells, lattice_constant, ...
                     thickness, epsilon_array, fill_factor] = wall_properties{:};
                        
                obj.add_grating_array(num_cells, lattice_constant, ...
                    thickness, epsilon_array, fill_factor, y_wall_center)  
            end
            
            %% extract the location of the core centers
            for i = 1:length(y_centers)-1
                ywall1 = obj.grating_properties{i}.y_pos;
                ywall2 = obj.grating_properties{i+1}.y_pos;
                gp = obj.grating_properties;
                %get mean 
                obj.core_centers_y(i) = (ywall1+ywall2)/2;
                obj.wall_bounds(i,:) = [gp{i}.wall_coords(2), gp{i+1}.wall_coords(1)];
            end
                
        end
        
      
   end
end