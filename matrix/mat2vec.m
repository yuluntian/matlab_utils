% function v = mat2vec(M)
% Vectorize a matrix by concatenating columns into a vector
%
% Yulun Tian
function v = mat2vec(M)
    v = reshape(M, [], 1);
end