% function T = extract_blktree_matrix(A, p)
% Given a p-by-p block-structured matrix symmetric matrix A, 
% returns a sparse matrix T that consists of the diagonal blocks of A, 
% together with off-digonal blocks that form a tree sparsity pattern. 
% The off-diagonal blocks are selected by computing the maximum spanning 
% tree where each edge is weighted by the norm of the corresponding block. 
%
% Yulun Tian
function X = extract_blktree_matrix(A, p)
n = size(A,1);
assert(n == size(A,2));
num_blocks = n / p;

% First, extract the diagonal blocks of A
X = extract_blkdiag_matrix(A, p);

% Next, select off-diagonal blocks
edge_src = [];
edge_dst = [];
edge_weight = [];
for i = 1:num_blocks
    for j = i+1:num_blocks
        idxs = (i-1)*p + 1 : i*p;
        jdxs = (j-1)*p + 1 : j*p;
        w = norm(A(idxs, jdxs), 'fro');
        if w > 1e-8
            edge_src(end+1) = i;
            edge_dst(end+1) = j;
            edge_weight(end+1) = 1/w;  % we use 1/w as the weight because our implementation computes the mininum spanning tree
        end
    end
end
G_ = graph(edge_src, edge_dst, edge_weight);
T_ = minspantree(G_, 'Type', 'forest');
[s, t] = findedge(T_);
for k = 1:length(s)
    i = s(k);
    j = t(k);
    idxs = (i-1)*p + 1 : i*p;
    jdxs = (j-1)*p + 1 : j*p;
    X(idxs, jdxs) = X(idxs, jdxs) + A(idxs, jdxs);
    X(jdxs, idxs) = X(jdxs, idxs) + A(jdxs, idxs);
end


end