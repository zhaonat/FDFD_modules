
function [inv_tensor] = inv_block(eps_tensor)
    
    %% calculate sub-blocks from inverting a 2x2 block matrix
    % hopefully the sub-blocks are all diagonal (which they should be)

    %% unravel
    Texx = eps_tensor{1,1}; %A
    Texy = eps_tensor{1,2}; %B
    Teyx = eps_tensor{2,1}; %C
    Teyy = eps_tensor{2,2}; %D
    
    inv_tensor_11 = Texx^-1 + Texx^-1*Texy*(Teyy-Teyx*Texx^-1*Texy)^-1*Teyx*Texx^-1;
    inv_tensor_12 = -Texx^-1*Texy*(Teyy-Teyx*Texx^-1*Texy)^-1;
    inv_tensor_21 = -Teyy^-1*Teyx*(Texx-Texy*Teyy^-1*Teyx)^-1;
    inv_tensor_22 = Teyy^-1+Teyy^-1*Teyx*(Texx-Texy*Teyy^-1*Teyx)^-1*Texy*Teyy^-1;

    inv_tensor = {inv_tensor_11, inv_tensor_12; inv_tensor_21, inv_tensor_22};
    
    
end