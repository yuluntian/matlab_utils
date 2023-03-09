% function [gi, gj, Hii, Hij, Hjj] = differentiate_relative_pose_measurement_2d( ...
%    Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)
% Compute the Riemannian gradient and optionally the Riemannian Hessian of
% a relative pose measurement between (Ri, ti) and (Rj, tj).
%
% Yulun Tian
function [gi, gj, Hii, Hij, Hjj] = differentiate_relative_pose_measurement_2d( ...
    Ri, ti, Rj, tj, Rij, tij, kappa, tau, options)

% Differentiate relative rotation measurement
if nargout == 2
    [gi_r, gj_r] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options);
    gi = [gi_r; zeros(2,1)];
    gj = [gj_r; zeros(2,1)];
else
    [gi_r, gj_r, Hii_r, Hjj_r, Hij_r] = differentiate_relative_rotation_measurement_2d(Ri, Rj, Rij, kappa, options);
    gi = [gi_r; zeros(2,1)];
    gj = [gj_r; zeros(2,1)];
    Hii = zeros(3,3); Hii(1,1) = Hii_r;
    Hjj = zeros(3,3); Hjj(1,1) = Hjj_r;
    Hij = zeros(3,3); Hij(1,1) = Hij_r;
end

G = gen2;
% Compute gradient of relative translation measurement
rgrad = [trace(G'*Ri'*(ti*tij' - tj*tij'));
         ti - tj + Ri * tij;
         0;
         tj - ti - Ri * tij];
rgrad = tau * rgrad;
gi = gi + rgrad(1:3);
gj = gj + rgrad(4:6);

% Compute Hessian of relative translation measurement, if requested
if nargout > 2
    A = trace(G'*G'*Ri'*(ti*tij' - tj * tij'));
    B = Ri * G * tij;
    %O = zeros(2,2);
    o = zeros(2,1);
    I = eye(2);
    rHess = [A  B' 0  -B';
             B  I  o  -I;
             0  o' 0   o';
            -B -I  o   I];
    rHess = tau * rHess;
    Hii = Hii + rHess(1:3,1:3);
    Hjj = Hjj + rHess(4:6,4:6);
    Hij = Hij + rHess(1:3,4:6);
end

end