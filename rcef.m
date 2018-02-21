function [U, monb] = rcef(A, tol)
% U - reduced column echelon form
% monb - monomial basis
%
% The reduction is performed by Gaussian elimination with
% column pivoting, based on Matlab's RREF routine

[m, n] = size(A);
% Compute the default tolerance if none was provided.
if (nargin < 2)
    tol = max(m,n)*eps(class(A))*norm(A,'inf');
end

% Loop over the entire matrix.
i = 1;
j = 1;
monb = [];
U = A;

while ((i <= m) && (j <= n))
    % Find value and index of largest element in the remainder of row j
    [p,k] = max(abs(U(i, j:n))); k = k+j-1;
    if (p <= tol)
        % The row is negligible, zero it out.
        U(i, j:n) = zeros(1, n-j+1);
        i = i + 1;
    else
        % Remember row index
        monb = [monb i];
        % Swap j-th and k-th columns
        U(i:m, [j k]) = U(i:m, [k j]);
        % Divide the pivot column by the pivot element
        U(i:m, j) = U(i:m, j) / U(i, j);
        % Subtract multiples of the pivot column from all the other columns
        for k = [1:j-1 j+1:n]
            U(i:m, k) = U(i:m, k) - U(i, k)*U(i:m, j);
        end
        i = i + 1;
        j = j + 1;
    end
    
end

end