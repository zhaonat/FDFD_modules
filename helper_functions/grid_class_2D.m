classdef grid_class_2D
   
   % attempt to standardize and abstract the grid
   % for use in simulations
    
   properties
      grid
      N
      Npml
      xrange
      yrange
      L
      Lpml
      shapes
      
   end
   
   methods
       
      function obj = grid_class_2D(eps, N, Npml, xrange, yrange, Lpml)
          obj.grid = eps;
          obj.N = N;
          obj.Npml = Npml;
          obj.xrange = xrange;
          obj.yrange = yrange;
          obj.L = [diff(xrange), diff(yrange)];
          obj.Lpml  = Lpml;    
      end

      function r = getCenter(obj)
         r = [round(obj.N(1)/2), round(obs.N(2)/2)];
      end

      function [nxpml, nypml] = getPML(obj)
         nxpml = [obj.Npml(1), obj.N(1)-obj.Npml(1)];
         nypml = [obj.Npml(2), obj.N(2)-obj.Npml(2)];
      end
      
      
   end
   
end