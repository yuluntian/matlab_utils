% function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options)
% Compute the Riemannian gradient and optionally the Riemannian Hessian of
% a relative measurement between two 2D rotations Ri and Rj.
% 
% 
% Yulun Tian
function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options)

% Compute error residual
RError = (Ri * Rij)' * Rj;
thetaij = atan2(RError(2,1), RError(1,1));

% Compute Riemannian gradient
grad_angular_dist = [-1; 1];   
if strcmp(options.rotation_distance, 'chordal')
    fdot_theta = 2 * sin(thetaij);
    fddot_theta = 2 * cos(thetaij);
else
    error('Unknown rotation distance: %s!', options.rotation_distance)
end
rgrad = fdot_theta * kappa * grad_angular_dist;
gi = rgrad(1);
gj = rgrad(2);

if nargout > 2
    % TODO: derive
    H12 = fddot_theta;
    Hii = H12 * kappa;
    Hjj = H12 * kappa;
    Hij = -H12 * kappa;
end

end