% function w = log2d(R)
% Return the scalar that corresponds to the element of the lie algebra
% so(2) for the input 2D rotation matrix R
%
% Yulun Tian
function w = log2d(R)
assert(size(R,1) == 2);
assert(size(R,2) == 2);
W = logm(R);
w = W(2,1);
end