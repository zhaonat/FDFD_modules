function [inv_cell] = inv_block_matrix(block_matrix)
%INV_BLOCK_MATRIX Summary of this function goes here
%   Detailed explanation goes here
% inverts a 2x2 block matrix where the input is a 2x2 cell

A = block_matrix{1,1};
B = block_matrix{1,2};
C = block_matrix{2,1};
D = block_matrix{2,2};

inv_cell{1,1} = inv(A-B*(D\C));
inv_cell{1,2} = -(A\B)\((D-C*(A\B)));
inv_cell{2,1} = -(D\C)\(A-B*(D\C));
inv_cell{2,2} = inv(D-C*(A\B));

end

