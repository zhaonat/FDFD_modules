classdef grid_class_2D
   
   % attempt to standardize and abstract the grid
   % for use in simulations
    
   properties
      grid % since epsilon is part of the properties, do we need to know Npml?
      N
      Npml
      xrange
      yrange
      L
      Lpml
      shapes %this is a cell which contains metadata on things added to the grid...
      
      % cell entry = cell with a {'shape name', metadata1, metadata2, etc}
      % doesn't seem like it's enough.. if I see a cell array in matlab...
      % I don't know shit unless I look up the original function..
      
   end
   
   methods
       
      function obj = grid_class_2D(N, Npml, xrange, yrange, Lpml)
          obj.grid = ones(N); %we shouldn't allow the user to put it in an epsilon to init...otherwise there may be shit on the grid
                              %without any metadata
          obj.N = N;
          obj.Npml = Npml;
          obj.xrange = xrange;
          obj.yrange = yrange;
          obj.L = [diff(xrange), diff(yrange)];
          obj.Lpml  = Lpml;   
          obj.shapes = cell(1);
      end

%       function r = getCenter(obj) %get center of the grid
%          r = [round(obj.N(1)/2), round(obs.N(2)/2)];
%       end

      function [nxpml, nypml] = getPML(obj)
         nxpml = [obj.Npml(1), obj.N(1)-obj.Npml(1)];
         nypml = [obj.Npml(2), obj.N(2)-obj.Npml(2)];
      end
      
      function[mask] = get_shape_mask(obj, background_epsilon)
          %not sure how useful this is.
          if(nargin <1)
             background_epsilon = 1; 
          end
          mask = obj.grid > background_epsilon; % in general, not enough to know what shapes...  
      end
      
      
   end
   
end