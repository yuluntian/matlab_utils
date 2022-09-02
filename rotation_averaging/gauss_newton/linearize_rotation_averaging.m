% function [residual, J, W] = linearize_rotation_averaging(measurements, R, options)
% Compute residual and Jacobian for 2D/3D rotation averaging problem.
function [residual, J, W] = linearize_rotation_averaging(measurements, R, options)

d = length(measurements.t{1});
n = max(max(measurements.edges));
m = size(measurements.edges, 1);
assert(strcmp(options.rotation_distance, 'chordal'));

if d == 2
    p = 1;  % state dimension
    rdim = 4;  % residual dimension
elseif d == 3
    p = 3; 
    rdim = 9;
else
    error('Invalid dimension %i', d);
end

residual = [];
Jrows = [];
Jcols = [];
Jvals = [];
Wvals = [];
for k = 1:m
    %% Read k-th measurement
    i = measurements.edges(k,1);
    j = measurements.edges(k,2);
    idxs = ((i-1)*d+1) : i*d;
    jdxs = ((j-1)*d+1) : j*d;
    Ri = R(:,idxs);
    Rj = R(:,jdxs);
    kappa = measurements.kappa{k};
    Rij = measurements.R{k};

    index_offset_res = (k-1) * rdim;
    index_offset_i = (i-1) * p;
    index_offset_j = (j-1) * p;

    % Linearize relative pose measurement
    [rvec, Ji, Jj] = linearize_relative_rotation_measurement(Ri, Rj, Rij, options);
    residual = [residual; rvec];

    % Push Ji
    for r = 1:rdim
        for c = 1:p
            Jrows(end+1) = index_offset_res + r;
            Jcols(end+1) = index_offset_i + c;
            Jvals(end+1) = Ji(r,c);
        end
    end

    % Push Jj
    for r = 1:rdim
        for c = 1:p
            Jrows(end+1) = index_offset_res + r;
            Jcols(end+1) = index_offset_j + c;
            Jvals(end+1) = Jj(r,c);
        end
    end

    % Push weights
    Wvals = [Wvals; kappa * ones(rdim, 1)];

end

% Assemble sparse Jacobian and weight matrices
J = sparse(Jrows, Jcols, Jvals, rdim * m, p * n);
W = spdiags(Wvals, 0, rdim * m, rdim * m);

end