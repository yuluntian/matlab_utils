% function W = hat3(w)
% Convert a 3-vector to the corresponding skew-symmetric matrix that
% represents the tangent vector of SO(3)
function W = hat3(w)

assert(size(w,1) == 3);
assert(size(w,2) == 1);
W = [0 -w(3) w(2);
     w(3) 0 -w(1);
     -w(2) w(1) 0];


end