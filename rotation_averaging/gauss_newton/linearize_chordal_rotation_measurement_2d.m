% function [r, Ji, Jj] = construct_relative_rotation_jacobian_2d(Ri, Rj, Rij, options)
% Compute the residual and Jacobian of a 2d relative rotation measurement
%
% Input:
% options.tangent_space_parametrization: 'global' or 'local'
%
% Output:
% r: residual vector is the vectorized version of (Rj - Ri*Rij)
% Ji: jacobian matrix of Ri
% Jj: jacobian matrix of Rj
%
% Yulun Tian
function [r, Ji, Jj] = linearize_chordal_rotation_measurement_2d(Ri, Rj, Rij, options)

% Compute residual
rMatrix = Rj - Ri * Rij;
r = reshape(rMatrix, [], 1);

% Compute Jacobian
G = [0 -1; 1 0];  % generator on so(2)
if strcmp(options.tangent_space_parametrization, 'global')
    Ji = reshape(- G * Ri * Rij, [], 1);
    Jj = reshape(G* Rj, [], 1);
elseif strcmp(options.tangent_space_parametrization, 'local')
    Ji = reshape(-Ri * G * Rij, [], 1);
    Jj = reshape(Rj * G, [], 1);
else
    error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
end

end