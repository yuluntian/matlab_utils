% function [g, H] = differentiate_pgo(measurements, R, t, options)
% Construct the Riemannian gradient and optionally the Riemannian Hessian
% of a 3D pose graph optimization problem.
%
% Yulun Tian
function [g, H] = differentiate_pgo(measurements, R, t, options)

if nargin < 4
    options = struct;
end
if ~isfield(options, 'rotation_distance')
    options.rotation_distance = 'chordal';
end
d = size(measurements.t{1},1); 
assert(d==2 || d == 3);
if d == 2
    diff_func = @differentiate_relative_pose_measurement_2d;
    p = 3;
else
    diff_func = @differentiate_relative_pose_measurement;
    p = 6;
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
    ti = t(:,i);
    tj = t(:,j);
    kappa = measurements.kappa{k};
    tau = measurements.tau{k};
    Rij = measurements.R{k};
    tij = measurements.t{k};
    
    % Add gradient and Hessian contributed by this measurement
    idxs_ = ((i-1)*p+1) : i*p;
    jdxs_ = ((j-1)*p+1) : j*p;
    if nargout == 1
        [gi, gj] = diff_func(Ri, ti, Rj, tj, Rij, tij, kappa, tau, options);
    else
        [gi, gj, Hii, Hij, Hjj] = diff_func(Ri, ti, Rj, tj, Rij, tij, kappa, tau, options);
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

% Construct sparse Hessian matrix
if nargout > 1
    H = sparse(Hrows, Hcols, Hvals, p * n, p * n);
    H = (H + H') / 2;
end


end