% function J = jacobian3(phi)
% compute the left Jacobian matrix 
% input phi is a 3-vector that corresponds to the element in so(3)
%
% Reference: Barfoot equation (7.35a)
%
% Yulun Tian
function J = jacobian3(phi)
assert(size(phi,1) == 3);
assert(size(phi,2) == 1);
angle = norm(phi);

if angle < 1e-12
    J = eye(3);
else
    axis = phi / angle;
    J = (sin(angle)/angle) * eye(3) + (1- sin(angle)/angle) * (axis *axis') + ((1-cos(angle))/angle) * hat3(axis);
end

end