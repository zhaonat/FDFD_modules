
function [filter] = create_filter(N, edge_width,alpha)
    %% Creates a filter to zero out the near-field at the edges
    %% for structures which are infinitely extended in one direction
    if(nargin<3)
       alpha = 0.5; 
    end
    filter = ones(1,N);
    
    x = linspace(-10,10, edge_width);
    y = (1+tanh(alpha*x))/2;
    filter(1:edge_width) = y;
    filter(end-edge_width+1:end) = fliplr(y);
    
    

end