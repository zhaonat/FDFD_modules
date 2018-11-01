function Dws = createDws_dirichlet_2D(w, s, dL,N)

    %% more transparent way to construct the dws matrices...but only for 2D
    %% MEANT FOR CREATING A PEC< BUT CURRENTLY DOESN'T WORK
    %% theoretically, doing a mask versus constructing the operator shouldn't be different
    
    Nx = N(1); Ny = N(2);
    Ix = speye(Nx); 
    Iy = speye(Ny); 
    
    %% this is crucial
    Ix(1,1) = 0; Ix(Nx,Nx) = 0;
    Iy(1,1) = 0; Iy(Ny,Ny) = 0;

    %% Create derivative operators
    switch w
        case 'x'
            if s == 'f'
                dxf = -Ix + circshift(Ix, [0 1]);
                %zero out anything connected to the edges...
                dxf(Nx,:) = 0; dxf(1,:) = 0;
                dxf(:,Nx) = 0; dxf(:,1) = 0;
                Dws = 1/dL(1) * kron(Iy, dxf); 
            else
                dxb = Ix - circshift(Ix, [0 -1]); 
                dxb(1,:) = 0;  dxb(Nx,:) = 0;
                dxb(:,1) = 0;   dxb(:,Nx) = 0;
                Dws = 1/dL(1) * kron(Iy, dxb); 
            end


        case 'y'
            if s == 'f'
                dyf = -Iy + circshift(Iy, [0 1]); 
                dyf(Ny,:) = 0; dyf(:,1) =0;  
                dyf(1,:) = 0; dyf(:,Ny) = 0;
                Dws = 1/dL(2) * kron(dyf, Ix); 
            else
                dyb = Iy - circshift(Iy, [0 -1]); 
                dyb(1,:) = 0;  dyb(:,Ny) = 0;         
                dyb(Ny,:) = 0; dyb(:,1) = 0;

                Dws = 1/dL(2) * kron(dyb, Ix); 
            end
    end   
end