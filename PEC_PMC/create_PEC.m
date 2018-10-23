function PEC_mask = create_PEC(N, w)

    %% HOW DO WE IMPLEMENT A PEC in the frequency domain
    % w = 'x', or 'y'
    % returns PEC_mask, which is an MxM operator which we can multiply with
    % the linear matrix A
    % all tangential electric fields go to 0
    % but normal electric fields are allowed to live (due to surface charge
    % induced)
    % this also means that normal H fields must disappear
    % the complexity is the offset fields in a Yee Cell
    % __Ex__
    %|      |  in the TM case, Ex or Ey may have to disappear
    %|  Hz  |  but tangential components are preserved
    %Ey     Ey 
    %|__Ex__|
    % __Hx__
    %|      |
    %|  Ez  |
    %Hy     Hy
    %|__Hx__|
    % ideally, we'd want control over which boundaries we specify the PEC
    % over: it appears, we must do left and right simultaneously
    % as in, we can't specify left wall to be PEC but the other to be
    % periodic
    % arguments: direction ('x' or 'y')
    % once we do this, it seems like the only difference is selective
    % changes in the derivative operators
    % i.e. a PEC, TM case at x = 0; Dxb becomes dirichlet
    % dirichlet, but Dyf, d
    M = prod(N);
    mask = ones(N);
    xn = 1:N(1);
    yn = 1:N(2);
    [Xn,Yn] = meshgrid(xn,yn);
    Xn = Xn.'; Yn = Yn.';
    if(w == 'x')
        mask(Xn == 1) = 0;
        mask(Xn == N(1)) =0;
    elseif(w == 'y')
        mask(Yn == 1) = 0;
        mask(Yn == N(2)) =0;   
    end
    % right now, if we wrap the entire grid in a PEC, the field pattern is
    % not symmetric...for sufficiently large domain size to wavelength
    % consequence of dispersion?
    
    PEC_mask = spdiags(mask(:), 0, M,M);

    
end
