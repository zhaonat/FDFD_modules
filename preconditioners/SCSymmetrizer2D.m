

function [Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxb, Syb)
    N = size(Sxf);

    sxf = diag(Sxf);
    syf = diag(Syf);

    numerator = sqrt(sxf).*sqrt(syf);
    denominator = 1./(numerator);
    Pr = spdiags(numerator, 0, N(1), N(2));
    Pl = spdiags(denominator, 0, N(1), N(2));
end