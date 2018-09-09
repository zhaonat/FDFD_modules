
function grading = logarithmic_grading(h0, hf,N)
    % alpha: grading factor
    % N, number of steps to grade down to
    grading = logspace(log10(h0), log10(hf), N);

end