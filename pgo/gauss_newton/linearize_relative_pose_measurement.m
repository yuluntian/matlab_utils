function [r, Ji, Jj] = linearize_relative_pose_measurement(Ri, ti, Rj, tj, Rij, tij, options)
d = size(Ri,1);
assert(d == 2 || d == 3);
assert(strcmp(options.rotation_distance, 'chordal'));

% Linearize the relative rotation measurement
[rvec_R, JR_Ri, JR_Rj] = linearize_relative_rotation_measurement(Ri,Rj,Rij,options);

% Compute translation residual
rvec_t = tj - ti - Ri * tij;

% Compute translation Jacobian
if d == 2
    Jt_tj = eye(2);
    Jt_ti = -eye(2);
    G = [0 -1; 1 0];  % generator on so(2)
    if strcmp(options.tangent_space_parametrization, 'global')
        Jt_Ri = - G * Ri * tij;
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Jt_Ri = - Ri * G * tij;
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end
else
    Jt_tj = eye(3);
    Jt_ti = -eye(3);
    if strcmp(options.tangent_space_parametrization, 'global')
        Jt_Ri = hat3(Ri * tij);
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Jt_Ri = Ri * hat3(tij);
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end

end

% Assemble overall residual and Jacobian
r = [rvec_R; rvec_t];
if d == 3
    Ji = [JR_Ri zeros(9,3);
          Jt_Ri Jt_ti];
    Jj = [JR_Rj      zeros(9,3);
          zeros(3,3) Jt_tj];
else
    Ji = [JR_Ri zeros(4,2);
          Jt_Ri Jt_ti];
    Jj = [JR_Rj      zeros(4,2);
          zeros(2,1) Jt_tj];
end

end