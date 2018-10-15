
function[] = band_plotter(k_vector_cell, omega_vector_cell)

    %% we want this function to encapsulate an idea about TRACKING BANDS
    %, which we realize with color coding...
    num_samples = length(k_vector_cell)
    
    %the issue is that plotting the bands in arbitrary order doesn't
    %necessarily mean we plot band 1's point first, band 2's point
    %2nd...they could be scrambled.
    
    %the ordering we follow will be dependent on the type of band-structure
    %we compute...whether k was the eigenvalue or omega.
    
    for i = 1:num_samples
        k = k_vector_cell{i};
        omega = omega_vector_cell{i};
        num_bands = length(k);
        
    end

end