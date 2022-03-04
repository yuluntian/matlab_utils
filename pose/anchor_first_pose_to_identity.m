% function [R, t] = anchor_first_pose_to_identity(R0, t0)
% helper function to apply a global transformation to the input set of
% poses, such that the resulting first pose is at the identity.
function [R, t] = anchor_first_pose_to_identity(R0, t0)

d = size(t0,1);
n = size(t0,2);
assert(size(R0,2) == d*n);
assert(size(R0,1) == d);
T_src_dst.R = R0(:, 1:d);
T_src_dst.t = t0(:, 1);
T_dst_src = invert_pose_struct(T_src_dst);

t = transform_translations(t0, T_dst_src);
R = transform_rotations(R0, T_dst_src);

end