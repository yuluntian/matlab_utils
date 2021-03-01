% function th = rotm2d_to_euler(R)
% Compute the Euler angle (radian) from a given 2D rotation matrix
%
% Yulun Tian
function th = rotm2d_to_euler(R)
check_rotation_matrix(R);
assert(size(R,1) == 2);
th = atan2(R(2,1), R(1,1));
end

