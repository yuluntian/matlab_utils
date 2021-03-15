% function R_dst = transform_rotations(R_src, T_dst_src)
% Transform a set of rotations from the input 'src' frame to the output 
% 'dst' frame by applying the known transformation.
%
% Each set of translation is a d-by-dn matrix R = [R1 ... Rn]
% 
% Yulun Tian 
function R_dst = transform_rotations(R_src, T_dst_src)

d = size(R_src, 1);
assert(size(T_dst_src.R, 1) == d);

R_dst = T_dst_src.R * R_src;

end