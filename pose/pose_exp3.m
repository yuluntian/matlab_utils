% function T = pose_exp3(x)
% Compute the exponential mapping for SE(3) 
% x is a 6D vector x = [rho phi],
% where rho corresponds to translation, and phi corresponds to rotation
%
% Yulun Tian
function T = pose_exp3(x)

assert(size(x,1) == 6);
assert(size(x,2) == 1);
rho = x(1:3);
phi = x(4:6);

C = exp3(phi);
J = jacobian3(phi);
T = [C               J * rho;
       zeros(1,3)  1];

end