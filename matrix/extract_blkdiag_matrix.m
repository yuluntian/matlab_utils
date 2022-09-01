% function [D, blocks] = extract_blkdiag_matrix(A, bdims)
% Given a block-structured square matrix A with block size specified by
% bdims(1), ..., bdims(n), return matrix D consisting of the diagonal blocks of A.
%
% Yulun Tian
function [D, blocks] = extract_blkdiag_matrix(A, bdims)
n = size(A,1);
assert(n == size(A,2));
assert(n == sum(bdims));
num_blocks = length(bdims);
blocks = cell(1, num_blocks);

for block_id = 1:num_blocks
    offset = sum(bdims(1:(block_id-1)));
    d = bdims(block_id);
    idxs = (offset + 1) : (offset + d);
    blocks{block_id} = A(idxs, idxs);
end

D = sparse(blkdiag(blocks{:}));

end