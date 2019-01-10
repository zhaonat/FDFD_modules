
%% function which sorts eigenvalues and eigenmodes by parity (even or odd)

function [even_k_eigs, odd_k_eigs] = ...
    parity_filter_1D(eigen_modes, eigenvals)
    % input is a cell, where each cell contains an eigenmode profile

    Nx = length(eigen_modes{1});
    even_k_eigs = []; omega_even_eigs = [];
    odd_k_eigs = [];  omega_odd_eigs = [];
    
%     odd_bandstructure = cell(1);
%     even_bandstructure = cell(1);
    
    for i = 1:length(eigen_modes)
        first_half = real(eigen_modes{i}(1:floor(Nx/2)));
        second_half = real(eigen_modes{i}(floor(Nx/2)+1:end));
        if((norm(abs(first_half-flipud(second_half)))<0.02))
            even_k_eigs= [even_k_eigs, eigenvals(i)];
        else
            odd_k_eigs= [odd_k_eigs, eigenvals(i)];
        end
    end
    
end