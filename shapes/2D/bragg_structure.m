
classdef bragg_structure <handle %need this so we can update by refernece
   
   % every structure class needs to know N, xrange, yrange from the
   % original grid
   
   % this formalism does make it harder to add in multiple different
   % structures. Do we assume that creating a bragg_grid obj 
   
   properties
      N
      xrange
      yrange
      epsilon % grid
      d1 %sublayer thickness 1
      d2 %sublayer thickness 2
      n_layers %each layer is two sub-layers thick
      x_pos %x1, x2 x boundaries
      y_pos % offset from the bottom of the ygrid
      
      %% defining parameters for any bragg grid on a rectangular domain
      
   end
   
   methods
        %two methods: one, initialize the object
        %             2) add the structure
        function [obj] = bragg_structure(xrange, yrange, N)
           obj.N = N;
           obj.xrange = xrange;
           obj.yrange=  yrange;
           obj.epsilon =ones(N);           
        end
        
        
        
        function [obj] = add_bragg(obj, rx, y_offset, n_layers, epsilon, d1, d2)
            %d1: thickness of layer 1 in microns
            %d2: thickness of layer 2 in microns
            % n_layers: each layer consists of d1 and d2 with e1 and e2
            % relative dielectrics
            
            % rx [x1, x2]; specifies start x and end x of the grid
            % delta_y;     offset from the bottom fo the domain
            
            obj.d1 = d1;
            obj.d2 = d2;
            obj.n_layers = n_layers;
            obj.x_pos = rx;
            obj.y_pos = y_offset;
            
            L = [diff(obj.xrange), diff(obj.yrange)];
            dL = L./obj.N;
            dy = dL(2);
            
            % convert contexts of rx into node_coords
            [x_start,~] = coord_to_grid([rx(1),0], obj.N, obj.xrange, obj.yrange);
            [x_end,~ ] = coord_to_grid([rx(2),0], obj.N, obj.xrange, obj.yrange);
            % get x_start and x_end
            ny_offset = round(y_offset/dy);

            % get thicknesses of each layer in the bragg stack in terms of grid
            % pointsd
            nd1 = round((d1)/dy);
            nd2 = round((d2)/dy);
            nd = [nd1,nd2];
            n_lay = nd1+nd2;
            
            % start constructing the layers
            for i = 1:n_layers
                y_start =ny_offset+(i-1)*n_lay+1;
                for j = 1:2
                    y_end = y_start+nd(j);
                    obj.epsilon(x_start:x_end, y_start:y_end) = epsilon(j);
                    y_start = y_end;
                end

            end

        end       
   end
end