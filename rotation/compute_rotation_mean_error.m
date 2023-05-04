% Given two sets of rotations, first align the rotations and then compute
% the average error in degree
%
% Each set of rotations is represented as a d-by-dn block-structured matrix
%
% Yulun Tian
function [error_degree, err_deg_vector]= compute_rotation_mean_error(Ropt, R)

d = size(R,1);
n = size(R,2) / d;

err_deg_vector = zeros(1,n);
[~, R_aligned] = align_rotations(R, Ropt);
for i = 1:n
    idxs = (i-1)*d+1 : i*d;
    R1 = R_aligned(:,idxs);
    R2 = Ropt(:, idxs);
    err_chord = norm(R1-R2, 'fro');
    err_deg_vector(i) = rad2deg(chordal_to_angular(err_chord));
end


error_degree = mean(err_deg_vector);

end