

function [adjoint_field] = adjoint_field(A,grad_u)
    
    %for optimizing the field at a single point, grad_u is exactly a point
    %source, equal to the original grid_field
    
    adjoint_field = A\grad_u;

end