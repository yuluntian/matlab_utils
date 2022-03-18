% function Ad = pose_adj3(T)
% Compute the adjoint of the input 3D pose T
% Input T is the 4-by-4 pose matrix
%
% Yulun Tian
function Ad = pose_adj3(T)

assert(size(T,1) == 4);
assert(size(T,2) == 4);

C = T(1:3,1:3);
r  = T(1:3,4);
Ad = [C hat3(r) * C;
         zeros(3,3)  C];

end