function [U, T] = orderschur(X)
% [U,T] = ORDERSCHUR(X) computes the real Schur decomposition X = U*T*U'
% of a matrix X with real eigenvalues sorted in increasing order along
% the diagonal of T
%
% Algorithm: perform unordered Schur decomposition with Matlab's SCHUR
% function, then order the eigenvalues with Givens rotations, see
% [Golub, Van Loan. Matrix computations. 1996]

[U, T] = schur(X);
U = real(U); T = real(T);
n = size(X,1);

order = 0;
while ~order % while the order is not correct
    order = 1;
    for k = 1:n-1
        if T(k, k) - T(k+1, k+1) > 0,
            order = 0; % diagonal elements to swap
            % Givens rotation
            [c,s] = givens(T(k, k+1),T(k+1, k+1)-T(k, k));
            T(k:k+1, k:n) = [c s; -s c]'*T(k:k+1, k:n);
            T(1:k+1, k:k+1) = T(1:k+1, k:k+1)*[c s; -s c];
            U(1:n, k:k+1) = U(1:n, k:k+1)*[c s; -s c];
        end
    end
end

end