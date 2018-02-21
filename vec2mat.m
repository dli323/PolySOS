function X = vec2mat(x, r)

if nargin < 2
    r = floor(sqrt(length(x)));
    if r^2 ~= length(x)
        error('Argument X has to be a square matrix!')
    end
end

X = reshape(x, r, r);
end