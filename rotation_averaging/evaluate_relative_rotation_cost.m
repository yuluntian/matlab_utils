% function cost = evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options)
% Compute the cost of a single relative rotation measurement:
% f(Ri, Rj) = 0.5 * kappa * dist(Ri Rij, Rj)^2
% 
% Yulun Tian
function cost = evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options)

if strcmp(options.rotation_distance, 'chordal')
    cost = 0.5 * kappa * (norm(Ri * Rij - Rj, 'fro').^2);
elseif strcmp(options.rotation_distance, 'geodesic')
    d = size(Ri, 1); 
    assert(d == 2 || d == 3);
    if d == 2
        cost = 0.5 * kappa * norm(log2d(Ri * Rij * Rj')).^2;
    else
        cost = 0.5 * kappa * norm(log3d(Ri * Rij * Rj')).^2;
    end
else
    error('Unknown rotation distance: %s!', options.rotation_distance)
end

end