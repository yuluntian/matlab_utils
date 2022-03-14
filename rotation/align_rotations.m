% function [R_dst_src, R_src_aligned] = align_rotations(R_src, R_dst)
% Given two sets of rotations expressed in src and dst frames, find the
% global rotation R_dst_src, which transforms rotations in the src frame to
% the dst frame and minimizes the total error in chordal distance.
% 
% Each set of rotations is specified as a d-by-dn matrix, e.g., R_src = 
% [R_src_1 ... R_src_n].
%
% Reference: Theorem 5 in SE-Sync paper
% https://arxiv.org/pdf/1612.07386.pdf
% 
% Author: Yulun Tian
function [R_dst_src, R_src_aligned] = align_rotations(R_src, R_dst)

d = size(R_src, 1);
n = size(R_src, 2) / d;
assert(size(R_dst, 1) == d);
assert(size(R_dst, 2) == d*n);

X = R_dst;
Y = R_src;
[U,~,V] = svd(X*Y');
E = eye(d);
E(d,d) = det(U*V');
R_dst_src = U * E * V';
R_src_aligned = R_dst_src * R_src;

end