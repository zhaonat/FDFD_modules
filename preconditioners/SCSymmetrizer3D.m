

function [Pl, Pr] = SCSymmetrizer3D(Sxf, Syf, Szf, Sxb, Syb, Szb)
    N = size(Sxf);
    %% stretch factor
    sxf = diag(Sxf);
    syf = diag(Syf);
    szf = diag(Szf);

    numerator = sqrt((sxf.*syf.*szf));
    denominator = 1./(numerator);
    Pr = spdiags(numerator, 0, N(1), N(2));
    Pl = spdiags(denominator, 0, N(1), N(2));

end