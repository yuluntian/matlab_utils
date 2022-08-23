% Given two sets of rotations, first align the rotations and then compute
% the RMSE rotation error in degree.
%
% Each set of rotations is represented as a d-by-dn block-structured matrix
%
% Yulun Tian
function error_degree = compute_rotation_RMSE(Ropt, R)

d = size(R,1);
n = size(R,2) / d;

squared_error = compute_rotation_orbit_distance(Ropt, R).^2;

% root mean squared error (Chordal)
error_chordal = sqrt(squared_error/n);

% convert chordal distance to degree
error_degree = rad2deg(chordal_to_angular(error_chordal));

end