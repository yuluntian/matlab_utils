% function cost = evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options)
% Compute the cost of a single relative rotation measurement
% 
% Yulun Tian
function cost = evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options)

if strcmp(options.rotation_distance, 'chordal')
    cost = 0.5 * kappa * (norm(Ri * Rij - Rj, 'fro').^2);
else
    error('Unknown rotation distance: %s!', options.rotation_distance)
end

end