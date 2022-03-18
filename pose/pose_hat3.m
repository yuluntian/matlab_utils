% function X = pose_hat3(x)
% Convert the 6D vector x into the corresponding element in se(3)
% x = [rho phi], where rho corresponds to translation and phi corresponds
% to rotation.
%
% Yulun Tian
function X = pose_hat3(x)
assert(size(x,1) == 6);
assert(size(x,2) == 1);
rho = x(1:3);
phi = x(4:6);
X = [hat3(phi) rho;
       zeros(1,3)  0];
end