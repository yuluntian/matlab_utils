% function [gi, gj, Hii, Hij, Hjj] = differentiate_relative_pose_measurement( ...
%    Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)
% Compute the Riemannian gradient and optionally the Riemannian Hessian of
% a relative pose measurement between (Ri, ti) and (Rj, tj).
%
% Yulun Tian
function [gi, gj, Hii, Hij, Hjj] = differentiate_relative_pose_measurement( ...
    Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)

% Differentiate relative rotation measurement
if nargout == 2
    [gi_r, gj_r] = differentiate_relative_rotation_measurement(Ri, Rj, Rij, kappa, options);
    gi = [gi_r; zeros(3,1)];
    gj = [gj_r; zeros(3,1)];
else
    [gi_r, gj_r, Hii_r, Hjj_r, Hij_r] = differentiate_relative_rotation_measurement(Ri, Rj, Rij, kappa, options);
    gi = [gi_r; zeros(3,1)];
    gj = [gj_r; zeros(3,1)];
    Hii = zeros(6,6); Hii(1:3,1:3) = Hii_r;
    Hjj = zeros(6,6); Hjj(1:3,1:3) = Hjj_r;
    Hij = zeros(6,6); Hij(1:3,1:3) = Hij_r;
end

% Compute gradient of relative translation measurement
rgrad = [hat3(tij)'*Ri'*(tj-ti);
         ti - tj + Ri * tij;
         zeros(3,1);
         tj - ti - Ri * tij];
rgrad = tau * rgrad;
gi = gi + rgrad(1:6);
gj = gj + rgrad(7:12);

% Compute Hessian of relative translation measurement, if requested
if nargout > 2
    A_ = hat3(Ri'*ti) * hat3(tij) - hat3(Ri'*tj)*hat3(tij);
    A = (A_ + A_')/2;
    B = Ri * hat3(tij);
    O = zeros(3,3);
    I = eye(3);
    rHess = [A -B' O   B';
            -B  I  O  -I;
             O  O  O   O;
             B -I  O   I];
    rHess = tau * rHess;
    Hii = Hii + rHess(1:6,1:6);
    Hjj = Hjj + rHess(7:12,7:12);
    Hij = Hij + rHess(1:6,7:12);
end

end