
function grading = logarithmic_grading(h0, hf,N)
    % alpha: grading factor
    % N, number of steps to grade down to
    % final grading is 2N+1 in shape
    % direction: 'u' or 'd'
    temp_ratio = (hf / h0)^(1/(N)); 
    grading = transpose(h0 * temp_ratio.^(0 : N-1)); 

end