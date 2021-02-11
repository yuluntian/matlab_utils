% Return the orbit distance between X, Y \in O(d)^n
% c.f. SE-Sync paper Appendix C.1
function distance = compute_orthogonal_orbit_distance(X, Y, d, n)
[~,S,~] = svd(X*Y');
nucnorm = sum(S(:));
distance = sqrt(2*d*n - 2*nucnorm);
end