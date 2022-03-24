% function [r, JRi, Jti, Jtj] = linearize_translation_measurement_2d(Ri, ti, tj, tij, options)
% Compute the residual and Jacobian of a 2d relative translation measurement
%
% r is the residual vector r = tj - ti - Ri * tij
function [r, JRi, Jti, Jtj] = linearize_translation_measurement_2d(Ri, ti, tj, tij, options)

% Compute residual
r = tj - ti - Ri * tij;

% Compute Jacobian
Jtj = eye(2);
Jti = -eye(2);
G = [0 -1; 1 0];  % generator on so(2)
if strcmp(options.tangent_space_parametrization, 'global')
    JRi = - G * Ri * tij;
elseif strcmp(options.tangent_space_parametrization, 'local')
    JRi = - Ri * G * tij;
else
    error('Unknown tangent space parametrization: %s', options.tangent_space_parametrization)
end

end