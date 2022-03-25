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
r = mat2vec(Rj - Ri * Rij);

% Compute Jacobian
if d == 2
    G = [0 -1; 1 0];  % generator on so(2)
    if strcmp(options.tangent_space_parametrization, 'global')
        Ji = mat2vec(- G * Ri * Rij);
        Jj = mat2vec(G* Rj);
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Ji = mat2vec(-Ri * G * Rij);
        Jj = mat2vec(Rj * G);
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end
else
    G1 = gen3(1);
    G2 = gen3(2);
    G3 = gen3(3);
    if strcmp(options.tangent_space_parametrization, 'global')
        Ji = -[mat2vec(G1*Ri*Rij) mat2vec(G2*Ri*Rij) mat2vec(G3*Ri*Rij)];
        Jj =  [mat2vec(G1*Rj) mat2vec(G2*Rj) mat2vec(G3*Rj)];
    elseif strcmp(options.tangent_space_parametrization, 'local')
        Ji = -[mat2vec(Ri*G1*Rij) mat2vec(Ri*G2*Rij) mat2vec(Ri*G3*Rij)];
        Jj =  [mat2vec(Rj*G1) mat2vec(Rj*G2) mat2vec(Rj*G3)];
    else
        error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
    end
end

end