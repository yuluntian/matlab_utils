% function [g, H] = differentiate_rotation_averaging(measurements, R, options)
% Construct the Riemannian gradient and optionally the Riemannian Hessian
% for a 3D rotation averaging problem
% 
% Yulun Tian
function [g, H] = differentiate_rotation_averaging(measurements, R, options, problem_data)

if nargin < 3
    options = struct;
end
if nargin < 4
    problem_data = struct;
end
if ~isfield(options, 'rotation_distance')
    options.rotation_distance = 'chordal';
end
if ~isfield(options, 'tangent_space_parametrization')
    options.tangent_space_parametrization = 'local';
end
assert(strcmp(options.tangent_space_parametrization, 'local') || ...
           strcmp(options.tangent_space_parametrization, 'global') );

% For chordal distance, if only gradient is needed and connection laplacian is also provided, 
% use the faster method
if nargout == 1 && strcmp(options.rotation_distance, 'chordal')
    if isfield(problem_data, 'ConLap')
        g = differentiate_rotation_averaging_chordal(measurements, R, options, problem_data);
        return
    else
        warning('For chordal distance, provide problem_data.ConLap to enable faster computation.')
    end
end

% Use the default method, which is capable of also computing the exact
% Hessian and works for geodesic distance. However, this implementation can
% be slow if the problem has many edges.
d = size(measurements.R{1}, 1);
assert(d == 2 || d == 3);
if d == 2
    diff_func = @differentiate_relative_rotation_measurement_2d;
    p = 1;
else
    diff_func = @differentiate_relative_rotation_measurement;
    p = 3;
end
n = max(max(measurements.edges));
m = size(measurements.edges, 1);
g = zeros(p*n, 1);
if nargout > 1
    % Initialize Riemannian Hessian
    Hrows = [];
    Hcols = [];
    Hvals = [];
end

for k = 1:m
    % Read k-th measurement
    i = measurements.edges(k,1);
    j = measurements.edges(k,2);
    idxs = ((i-1)*d+1) : i*d;
    jdxs = ((j-1)*d+1) : j*d;
    Ri = R(:,idxs);
    Rj = R(:,jdxs);
    kappa = measurements.kappa{k};
    Rij = measurements.R{k};
    
    % Add gradient and Hessian contributed by this measurement
    idxs_ = ((i-1)*p+1) : i*p;
    jdxs_ = ((j-1)*p+1) : j*p;
    if nargout == 1
        [gi, gj] = diff_func(Ri, Rj, Rij, kappa, options);
    else
        [gi, gj, Hii, Hjj, Hij] = diff_func(Ri, Rj, Rij, kappa, options);
        % Assemble Hii
        row_offset = (i-1)*p;
        col_offset = (i-1)*p;
        for r = 1:p
            for c = 1:p
                Hrows(end+1) = row_offset + r;
                Hcols(end+1) = col_offset + c;
                Hvals(end+1) = Hii(r,c);
            end
        end
        % Assemble Hjj
        row_offset = (j-1)*p;
        col_offset = (j-1)*p;
        for r = 1:p
            for c = 1:p
                Hrows(end+1) = row_offset + r;
                Hcols(end+1) = col_offset + c;
                Hvals(end+1) = Hjj(r,c);
            end
        end
        % Assemble Hij
        row_offset = (i-1)*p;
        col_offset = (j-1)*p;
        for r = 1:p
            for c = 1:p
                Hrows(end+1) = row_offset + r;
                Hcols(end+1) = col_offset + c;
                Hvals(end+1) = Hij(r,c);
            end
        end
        % Assemble Hji
        Hji = Hij';
        row_offset = (j-1)*p;
        col_offset = (i-1)*p;
        for r = 1:p
            for c = 1:p
                Hrows(end+1) = row_offset + r;
                Hcols(end+1) = col_offset + c;
                Hvals(end+1) = Hji(r,c);
            end
        end
    end
    % Assemble gradient
    g(idxs_) = g(idxs_) + gi;
    g(jdxs_) = g(jdxs_) + gj;
end

% Convert gradient to global tangent space (Lie algebra) if requested
if strcmp(options.tangent_space_parametrization, 'global')
    if d == 3
        % For SO(3), the global tangent vector eta is related to the local tangent vector w via a linear map (Adjoint): 
        % w = A * eta, where A is a block-diagonal matrix, with Aii = Ri^T
        A = sparse(d * n, d * n);
        for i = 1:n
            idxs = (i-1)*d+1 : i*d;
            A(idxs, idxs) = R(:, idxs)';
        end
    else
        % For SO(2), the linear map (Adjoint) is simply Identity, because
        % rotations commute in 2D.
        A = speye(n);
    end
    % The global tangent space gradient is computed via chain rule
    g = A' * g;
end

% Construct sparse Hessian matrix
if nargout > 1
    H = sparse(Hrows, Hcols, Hvals, p * n, p * n);
    if strcmp(options.tangent_space_parametrization, 'global')
        H = A' * H * A;
        [H, sym_error] = make_symmetric(H);
    else
        [H, sym_error] = make_symmetric(H);
    end
    assert(sym_error < 1e-8, ...
        'Hessian is not symmetric! Error: %f', sym_error);
end

end