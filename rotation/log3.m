% function w = log3(R)
% Compute the logarithmic mapping for a 3d rotation R
% Returns the 3-vector that corresponds to the element of so(3)
% with rotation angle less than pi
% 
% Reference: Sec 7.1, Barfoot state estimation for robotics
%
% Yulun Tian
function w = log3(R)
assert(size(R,1) == 3);
assert(size(R,2) == 3);

% find the eigenvector corresponding to eigenvalue 1
[axis, eigval] = eigs(R - eye(3), 1, 'smallestabs');
assert(abs(eigval) < 1e-10);

angle = acos((trace(R)-1)/2);
assert(angle < pi + 1e-10);
w = angle * axis;

% flip axis if needed
R_ = exp3(w);
if norm(R-R_,'fro') > 1e-8
    w = -w;
end
end