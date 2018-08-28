
function grading = logarithmic_grading(N, alpha, direction)
    % alpha: grading factor
    % N, number of steps to grade down to
    % final grading is 2N+1 in shape
    % direction: 'u' or 'd'
    switch direction
        case 'u'
            powers = 0:N;

        case 'd'
            powers = N:-1:0;
            
    end
    grading = (alpha.^powers).';
end