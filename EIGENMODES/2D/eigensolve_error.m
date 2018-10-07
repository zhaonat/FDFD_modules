

function error = eigensolve_error(A,eigenvector, eigenvalue)
    %determine how close an eigenvector, eigenvalue solution is to
    %minimizing the eigenvalue equation.
    error = norm(full(A*eigenvector-eigenvalue*eigenvector));
end