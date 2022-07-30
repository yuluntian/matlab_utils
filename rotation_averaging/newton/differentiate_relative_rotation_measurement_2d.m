% function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options)
% Compute the Riemannian gradient and optionally the Riemannian Hessian of
% a relative measurement between two 2D rotations Ri and Rj.
% f(R_i, R_j) = 0.5 * kappa * dist(Ri Rij, Rj)^2
% 
% 
% Yulun Tian
function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options)

% Compute error residual
RError = (Ri * Rij)' * Rj;
thetaij = log2d(RError);

% Compute Riemannian gradient
grad_angular_dist = [-1; 1];  % gradient of the 2D angular distance is constant
if strcmp(options.rotation_distance, 'chordal')
    fdot_theta = 2 * sin(thetaij);
    fddot_theta = 2 * cos(thetaij);
elseif strcmp(options.rotation_distance, 'geodesic')
    fdot_theta = thetaij;
    fddot_theta = 1;
else
    error('Unknown rotation distance: %s!', options.rotation_distance)
end
rgrad = fdot_theta * kappa * grad_angular_dist;
gi = rgrad(1);
gj = rgrad(2);

if nargout > 2
    % Ref: Tron (E.40) which is a result of chain rule (C.8)
    % Note that for SO(2), DLog(R)/||Log(R)|| is zero, and thus the second
    % term in E.40 disappears
    H12 = fddot_theta;
    Hii = H12 * kappa;
    Hjj = H12 * kappa;
    Hij = -H12 * kappa;
end

end