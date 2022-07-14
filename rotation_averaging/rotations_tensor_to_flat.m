% function Rflat = rotations_tensor_to_flat(R)
% Flatten a 3D tensor of multiple rotations (i.e., each R(:,:,i) is a rotation)
% to a 2D matrix Rflat = [R(:,:,1) ... R(:,:,n)]
%
% Author: Yulun Tian
function Rflat = rotations_tensor_to_flat(R)
Rflat = reshape(R, size(R,1), []);
end