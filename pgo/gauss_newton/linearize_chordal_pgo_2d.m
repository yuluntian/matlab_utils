% function [residual, J, W] = linearize_chordal_pgo_2d(measurements, R, t, options)
% Compute the residual and Jacobian of 2D PGO problem under the chordal
% distance.
%
% Yulun Tian
function [residual, J, W] = linearize_chordal_pgo_2d(measurements, R, t, options)

d = length(measurements.t{1});
n = max(max(measurements.edges));
m = size(measurements.edges, 1);
assert(d == 2);

residual = zeros(6*m, 1);
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
    ti = t(:,i);
    tj = t(:,j);
    kappa = measurements.kappa{k};
    tau = measurements.tau{k};
    Rij = measurements.R{k};
    tij = measurements.t{k};

    index_offset_Rij = (k-1) * 6;
    index_offset_tij = index_offset_Rij + 4;
    index_offset_Ri = (i-1) * 3;
    index_offset_ti = index_offset_Ri + 1;
    index_offset_Rj = (j-1) * 3;
    index_offset_tj = index_offset_Rj + 1;

    %% Linearize rotation measurement
    [res_R, Jk_i, Jk_j] = linearize_chordal_rotation_measurement_2d(Ri, Rj, Rij, options);
    % Push residual vector
    residual(index_offset_Rij+1 : index_offset_Rij+4) = res_R;
    % Push Jacobians
    for l = 1:4
        Jrows(end + 1) = index_offset_Rij + l;
        Jcols(end+1) = index_offset_Ri + 1;
        Jvals(end+1) = Jk_i(l);
    end
    for l = 1:4
        Jrows(end + 1) = index_offset_Rij + l;
        Jcols(end+1) = index_offset_Rj + 1;
        Jvals(end+1) = Jk_j(l);
    end
    % Push weights
   Wvals = [Wvals; kappa * ones(4,1)];

    %% Linearize translation measurement
    [res_t, JRi, Jti, Jtj] = linearize_translation_measurement_2d(Ri, ti, tj, tij, options);
    % Push residual vector
    residual(index_offset_tij+1: index_offset_tij + 2) = res_t;
    % Push Jacobians
    for l = 1:2
        Jrows(end+1) = index_offset_tij + l;
        Jcols(end+1) = index_offset_Ri + 1;
        Jvals(end+1) = JRi(l);
    end
    for r = 1:2
        for c = 1:2
            Jrows(end+1) = index_offset_tij + r;
            Jcols(end+1) = index_offset_ti + c;
            Jvals(end+1) = Jti(r,c);
        end
    end
    for r = 1:2
        for c = 1:2
            Jrows(end+1) = index_offset_tij + r;
            Jcols(end+1) = index_offset_tj + c;
            Jvals(end+1) = Jtj(r,c);
        end
    end
    % Push weights
    Wvals = [Wvals; tau * ones(2,1)];
end

% Assemble sparse Jacobian and weight matrices
J = sparse(Jrows, Jcols, Jvals, 6 * m, 3 * n);
W = spdiags(Wvals, 0, 6 * m, 6 * m);
end