% Return the orbit distance between X, Y \in SO(d)^n
% c.f. SE-Sync paper Appendix C.1
function distance = compute_rotation_orbit_distance(X, Y)
d = size(X,1);
n = size(X,2)/d;
assert(size(Y,1) == d);
assert(size(Y,2) == d*n);
[U,S,V] = svd(X*Y');
E = eye(d);
E(d,d) = det(U*V');
sq_dist = max(2*d*n - 2*trace(E*S), 0);
distance = sqrt(sq_dist);
end