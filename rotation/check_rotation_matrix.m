function check_rotation_matrix(R)
d = size(R,1);
assert(size(R,2) == d);

assert(norm(R'*R - eye(d), 'fro') < 1e-6);
assert(abs(det(R)-1) < 1e-6);
end