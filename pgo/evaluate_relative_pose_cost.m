% function cost = evaluate_relative_pose_cost(Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)
% Compute the cost of a single relative pose measurement
% 
% Yulun Tian
function cost = evaluate_relative_pose_cost(Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)

rotation_cost = evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options);
translation_cost = 0.5 * tau * (norm(tj - ti - Ri*tij).^2);
cost = rotation_cost + translation_cost;

end