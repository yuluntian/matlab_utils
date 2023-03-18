% function g = differentiate_rotation_averaging_chordal(measurements, R, options, problem_data)
% Construct the Riemannian gradient for a rotation averaging problem under
% the squared chordal distance
% 
% Yulun Tian
function g = differentiate_rotation_averaging_chordal(measurements, R, options, problem_data)
assert(strcmp(options.rotation_distance, 'chordal'));
assert(strcmp(options.tangent_space_parametrization, 'local') || ...
           strcmp(options.tangent_space_parametrization, 'global') );
assert(isfield(problem_data, 'ConLap'));
d = size(measurements.R{1}, 1);
assert(d == 2 || d == 3);
n = max(max(measurements.edges));
p = d * (d-1) / 2;
g = zeros(p*n, 1);   % Final gradient vector to be returned

% Obtain the manopt manifold object
if ~isfield(problem_data, 'Manifold')
    M = stiefelstackedfactory(n, d, d);
else
    M = problem_data.Manifold;
end

% Compute Riemannian gradient as the projection of Euclidean gradient to
% the tangent space.
euc_grad_T = problem_data.ConLap * R';
grad_T = M.proj(R', euc_grad_T);
G = grad_T';

% Convert gradient tangent vector
% from matrix form Gi = RiSi to vector form (in the standard basis)
for i = 1:n
    idxs = (i-1)*d + 1 : i*d;
    jdxs = (i-1)*p + 1 : i*p;
    Ri = R(:, idxs);
    Gi = G(:, idxs);
    % Si is the skew symmetric matrix that represents the tangent vector
    Si = Ri' * Gi;
    ss_err = norm(Si + Si', 'fro');
    if ss_err > 1e-6
        warning('Input is not skew-symmetric!');
    end
    if d == 3
        g(jdxs) = 2 * vee3(Si);  % the extra factor of 2 is due to the difference in the definition of inner product
    else
        g(jdxs) = 2* vee2(Si);
    end
    
    % Convert gradient to global parametrization, if requested
    if d == 3 && strcmp(options.tangent_space_parametrization, 'global') 
        g(jdxs) = Ri * g(jdxs);
    end
end





end