% function w = log3d(R)
% Compute the logarithmic mapping for a 3d rotation R
% Returns the 3-vector that corresponds to the element of so(3)
% with rotation angle less than pi
% 
% Reference: Ethan Eade, "Lie Groups for 2D and 3D Transformations"
%
% Yulun Tian
function w = log3d(R)
assert(size(R,1) == 3);
assert(size(R,2) == 3);

theta = acos((trace(R)-1)/2);
if theta < 1e-12
    w = zeros(3,1);
else
    w = vee3(theta * (R - R') / (2 * sin(theta)));
end


end