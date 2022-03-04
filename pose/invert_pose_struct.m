% function Tinv = invert_pose_struct(T)
% given a pose represented as a struct T, where T.R and T.t contains the
% rotation and translation vector, return a struct Tinv that represents the
% inverse of this pose.
%
% Author: Yulun Tian
function Tinv = invert_pose_struct(T)
R = T.R;
t = T.t;
d = size(t,1);
% Sanity check
assert(size(R, 1) == d);
assert(size(R, 2) == d);
assert(size(t, 2) == 1);
check_rotation_matrix(R);

% Tmat = [R,         t;
%         zeros(1,d) 1];
% Tinvmat = inv(Tmat);

Tinv = struct;
Tinv.R = R';
Tinv.t = -R' * t;

end