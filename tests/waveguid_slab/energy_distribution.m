%% function which extracts confined modes of a structure

% input is the set of vectors from the eigs function
%structure_mask is a binary map of 0's and 1's..1's indicate structure...
%U is reshaped
function ratio = energy_distribution(U, structure_mask)
    InsideStructure = U(structure_mask == 1);
    OutsideStructure = U(structure_mask == 0);
    
    energy_inside = sum(abs(InsideStructure).^2);
    energy_outside = sum(abs(OutsideStructure).^2);
    energy_total = energy_inside +energy_outside;
    ratio = energy_inside/(energy_total);
    
   
end