% function R = rotations_flat_to_tensor(Rflat)
% Convert a 2D matrix of multiple rotations
% Rflat = [R1 ... Rn] to a 3D tensor R,
% where R(:,:,i) = Ri.
%
% Author: Yulun Tian
function R = rotations_flat_to_tensor(Rflat)
R = reshape(Rflat, size(Rflat,1), size(Rflat,1), []);
end