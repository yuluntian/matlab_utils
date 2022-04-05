% function K = commutation_matrix(m,n)
% Return the commutation matrix for a matrix A of size m by n, such that 
% K vec(A) = vec(A')
% Ref: https://en.wikipedia.org/wiki/Commutation_matrix
%
% Yulun Tian
function K = commutation_matrix(m,n)

%[m, n] = size(A);

I = reshape(1:m*n, [m, n]);
I = I';
I = I(:);
K = eye(m*n);
K = K(I,:);

end