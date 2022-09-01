% function [D, blocks] = extract_blkdiag_matrix(A, block_size)
% Given a block-structured square matrix A, return matrix D consisting of the diagonal blocks of A.
% Input:
%   A: input matrix
%   block_size: if a scalar, it is assumed that all blocks have the same
%   size. If a vector, then block_size(i) is the dimension of the i-th
%   block.
%
% Yulun Tian
function [D, blocks] = extract_blkdiag_matrix(A, block_size)
n = size(A,1);
assert(n == size(A,2));
if length(block_size) > 1
    num_blocks = length(block_size);
    bdims = block_size;
else
    num_blocks = n / block_size;
    bdims = block_size * ones(1, num_blocks);
end
assert(n == sum(bdims));
blocks = cell(1, num_blocks);

for block_id = 1:num_blocks
    offset = sum(bdims(1:(block_id-1)));
    d = bdims(block_id);
    idxs = (offset + 1) : (offset + d);
    blocks{block_id} = A(idxs, idxs);
end

D = sparse(blkdiag(blocks{:}));

end