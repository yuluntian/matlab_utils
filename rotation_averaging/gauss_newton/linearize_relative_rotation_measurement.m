% function [r, Ji, Jj] = linearize_relative_rotation_measurement(Ri, Rj, Rij, options)
% Compute the residual and Jacobian of a relative rotation measurement
% between two 2D or 3D rotations.
%
% Input:
% options.rotation_distance: 'chordal' 
% options.tangent_space_parametrization: 'global' or 'local'
%
% Output:
% r: residual vector is the vectorized version of (Rj - Ri*Rij)
% Ji: jacobian matrix of Ri
% Jj: jacobian matrix of Rj
%
% TODO: support geodesic distance
%
% Yulun Tian
function [r, Ji, Jj] = linearize_relative_rotation_measurement(Ri, Rj, Rij, options)

d = size(Ri,1);
assert(d == 2 || d == 3);
assert(strcmp(options.rotation_distance, 'chordal'));

% Compute residual
rMatrix = Rj - Ri * Rij;
r = reshape(rMatrix, [], 1);

% Compute Jacobian
if d == 2
    G = [0 -1; 1 0];  % generator on so(2)
    if strcmp(options.tangent_space_parametrization, 'global')
        Ji = vecmat(- G * Ri * Rij);
        Jj = vecmat(G* Rj);
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Ji = vecmat(-Ri * G * Rij);
        Jj = vecmat(Rj * G);
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end
else
    G1 = gen3(1);
    G2 = gen3(2);
    G3 = gen3(3);
    if strcmp(options.tangent_space_parametrization, 'global')
        error('global parametrization for SO(3) is not implemented.');
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Ji = -[vecmat(Ri*G1*Rij) vecmat(Ri*G2*Rij) vecmat(Ri*G3*Rij)];
        Jj =  [vecmat(Rj*G1) vecmat(Rj*G2) vecmat(Rj*G3)];
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end
end

end