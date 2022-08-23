% function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement(Ri, Rj, Rij, kappa, options)
% Compute the Riemannian gradient and optionally the Riemannian Hessian of
% a relative measurement between two 3D rotations Ri and Rj:
% f(Ri, Rj) = 0.5 * kappa* dist(Ri Rij, Rj)^2
% 
% Note that the computed gradient and Hessian are in the local tangent
% space.
% 
% Yulun Tian
function [gi, gj, Hii, Hjj, Hij] = differentiate_relative_rotation_measurement(Ri, Rj, Rij, kappa, options)

% Compute error residual
RError = (Ri * Rij)' * Rj;
uij = log3d(RError);
thetaij = norm(uij);
if abs(thetaij) < 1e-12
    % Treat this error as zero
    thetaij = 0;
    uij = zeros(3,1);
    uij(1) = 1;
else
    uij = uij / thetaij;
end
O = zeros(3,3);
I = eye(3);
L = [Rij O;
     O   I];

% Compute Riemannian gradient
grad_angular_dist = [-uij; uij];   % (E.31) of Tron thesis
if strcmp(options.rotation_distance, 'chordal')
    fdot_theta = 2 * sin(thetaij);
elseif strcmp(options.rotation_distance, 'geodesic')
    fdot_theta = thetaij;
else
    error('Unknown rotation distance: %s!', options.rotation_distance)
end
rgrad = fdot_theta * kappa * L * grad_angular_dist;
gi = rgrad(1:3);
gj = rgrad(4:6);

% Compute Riemannian Hessian if requested
if nargout > 2
    alpha = alpha_func(thetaij, options);
    gamma = gamma_func(thetaij, options);
    beta  = beta_func(thetaij, options);
    Sij = alpha * eye(3) + gamma * (uij * uij');
    Aij = beta * hat3(uij);
    Hij_ = [Sij         (-Sij-Aij);
           (-Sij-Aij')  Sij];
    HR = kappa * L * Hij_ * L';
    HR = make_symmetric(HR);  % ensure Hessian is symmetric
    Hii = HR(1:3,1:3);
    Hjj = HR(4:6,4:6);
    Hij = HR(1:3,4:6);
end

end