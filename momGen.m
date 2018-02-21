function M = momGen(A, y)
% Generate the moment matrix from dual solution
% A - selection matrices in sedumi format
% y - dual solution of the SDP

v = -A'*y;
v(1,1) = 1; % Reconstruct the constant term
M = vec2mat(v);
end