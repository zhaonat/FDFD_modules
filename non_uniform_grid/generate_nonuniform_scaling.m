
function [dx_scale, dy_scale] = generate_nonuniform_scaling(Nft, drt)
    
    %Nft: 1st column is x, 2nd column is y
    %drt: list of discretizations...normalized by some reference
    % we can express drt as proportions of the largest discretization
    % available on the grid...but seems inefficient
    % advantage is we don't have to rewrite the pml sfactor

    Nx = sum(Nft(:,1));
    Ny = sum(Nft(:,2));
    dx_scale = ones(1,Nx);
    dy_scale = ones(1,Ny);
    
    num_regions = length(Nft(:,1));
    x0 = 1; y0 = 1;
    
    % Here, we can assume that all odd indices are fixed regions
    % even indices are transition regions
    
    for i = 1:2:num_regions
       dx_scale(x0:x0+Nft(i,1)-1) = drt(i,1);
       dy_scale(y0:y0+Nft(i,2)-1) = drt(i,2);
    

       if(i==num_regions) %no transition after last region
           x0 = x0+Nft(i,1);
           y0 = y0+Nft(i,1);
       else
           x0 = x0+Nft(i,1)+Nft(i+1,1);
           y0 = y0+Nft(i,1)+Nft(i+1,2);
       end
        
    end
    
    % do some sort of grading from region i to region i+1
    x0 = Nft(1,1); y0 = Nft(1,2);
    for i = 2:2:num_regions
        dx1 = drt(i-1,1); dx2 = drt(i+1,1);
        dy1 = drt(i-1,2); dy2 = drt(i+1,2);
        nxt = Nft(i,1); nyt = Nft(i,2);
        
        %need a function to grade smoothly from dr1 to dr2
        %equation to solve is there is some ration dr1/dr2 
        % need to multiply dr1 by constant nt times.

        grading_x = logspace(log10(dx1), log10(dx2), nxt+1);
        grading_y = logspace(log10(dy1), log10(dy2), nyt+1);

        dx_scale(x0:x0+nxt) = grading_x;
        dy_scale(y0:y0+nxt) = grading_y;
        x0 = x0+Nft(i,1)+Nft(i+1,1); 
        y0 = y0+Nft(i,2)+Nft(i+1,1);

    end
    

end