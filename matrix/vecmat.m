% function v = vecmat(M)
% Vectorize a matrix by concatenating columns into a vector
%
% Yulun Tian
function v = vecmat(M)
    v = reshape(M, [], 1);
end