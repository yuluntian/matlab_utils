% function [new_measurements, gt_info] = regenerate_measurements(measurements, xtrue, options)
% Given a measurement struct in SE-Sync format, regenerate each relative
% edge based on the specified noise model.
% 
% 'measurements' consists of:
% 
% edges:  An (mx2)-dimensional matrix encoding the edges in the measurement 
%      network; edges(k, :) = [i,j] means that the kth measurement is of the
%      relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%      that the states x_i are numbered sequentially as x_1, ... x_n.
% R:  An m-dimensional cell array whose kth element is the rotational part
%      of the kth measurement
% t:  An m-dimensional cell array whose kth element is the translational
%      part of the kth measurement
% kappa:  An m-dimensional cell array whose kth element gives the 
%      precision of the rotational part of the kth measurement. 
% tau:  An m-dimensional cell array whose kth element gives the precision
%      of the translational part of the kth measurement.
%
% 'xtrue' is a struct containing the 'ground truth' poses
% ---R: d-by-dn block-structured matrix of rotation estimates, i.e., R =
% [R1 R2 ... Rn]
% ---t: d-by-n block-structured matrix of translation vectors, i.e., t =
% [t1 t2 ... tn]
function [new_measurements, gt_info] = regenerate_measurements(measurements, xtrue, options)
d = length(measurements.t{1});
n = max(max(measurements.edges));  % number of poses
assert(d == 2 || d == 3, 'Invalid dimension: %g', d);
assert(size(xtrue.R, 1) == d);
assert(size(xtrue.R, 2) == d * n);
assert(size(xtrue.t, 1) == d);
assert(size(xtrue.t, 2) == n);

%% Parse options
rad_stddev = deg2rad(options.deg_stddev); % standard deviation of rotation measurements in radian
t_stddev = options.t_stddev;              % standard deviation of translational measurements in m
tau = 1/(t_stddev.^2);                    % translational measurement weights
kappa = (1/(rad_stddev.^2))/2;            % rotational measurement weights

% SUPPORT FOR SIMULATING OUTLIER MEASUREMENTS
% Percentage of outliers within all loop closures
if isfield(options, 'lc_outlier_prob')
    lc_outlier_prob = options.lc_outlier_prob;
    assert(lc_outlier_prob >= 0 && lc_outlier_prob <= 1);
else
    lc_outlier_prob = 0;
end

if isfield(options, 'lc_outlier_mode')
    lc_outlier_mode = options.lc_outlier_mode;
else
    lc_outlier_mode = 'replace';
end
assert(strcmp(lc_outlier_mode, 'add') || strcmp(lc_outlier_mode, 'replace'), ...
       'invalid outlier loop closure mode: s%. Options are add or replace.', ...
       lc_outlier_mode);
   
if isfield(options, 'lc_outlier_max_translation')
    randtmax = options.lc_outlier_max_translation;
else
    randtmax = 10;
end

%% Re-generate existing measurements
inlier_loop_closures = [];
outlier_loop_closures = [];

edges = measurements.edges;

for k = 1:size(edges,1)
    s = edges(k,1); % source index
    t = edges(k,2); % dest index
    
    Rsrc = xtrue.R(:, (s-1)*d+1 : s*d);
    Rdst = xtrue.R(:, (t-1)*d+1 : t*d);
    tsrc = xtrue.t(:, s);
    tdst = xtrue.t(:, t);
    
    % get ground truth transformation
    Ts = [Rsrc        tsrc;
          zeros(1,d)     1];
    Tt = [Rdst        tdst;
          zeros(1,d)     1];
    dT = inv(Ts)*Tt;
    dR = dT(1:d,1:d);
    dt = dT(1:d,end);

    % simulate noisy measurements
    if s + 1 == t
        % odometry edges are always inlier
        is_outlier = false;
    else
        if strcmp(lc_outlier_mode, 'replace')
            % replace existing loop closure with an outlier edge
            is_outlier = randsample([true false],1,true,...
                                    [lc_outlier_prob, 1-lc_outlier_prob]);
        else
            is_outlier = false;
        end
    end
    
    if ~is_outlier
        % Gaussian translation noise
        error_t = normrnd(0, t_stddev, [d,1]);
        
        % Langevin rotation noise
        if d == 3
            error_R = Langevin_sampler_3D(eye(3), kappa);
        else
            error_R = Langevin_sampler_2D(eye(2), kappa);
        end
        
        if s + 1 ~= t
            inlier_loop_closures = [inlier_loop_closures; s t];
        end
    else
        % Uniformly random translation noise 
        error_t = -randtmax + 2*randtmax*rand(d,1);
        
        % Uniformly random rotation noise
        if d == 3
            error_R = randrot(3,1);
        else
            th = pi * rand;
            error_R = [cos(th), -sin(th); 
                       sin(th), cos(th)];
        end
        outlier_loop_closures = [outlier_loop_closures; s t];
    end
    
    if ~is_outlier
        cell_t{k} = dt + error_t;
        cell_R{k} = dR * error_R;
    else
        cell_t{k} = error_t;
        cell_R{k} = error_R;
    end
    cell_tau{k} = tau;
    cell_kappa{k} = kappa;
end

%% In add mode, generate additional outlier loop closures
if strcmp(lc_outlier_mode, 'add')
    assert(isempty(outlier_loop_closures))
    
    num_outliers = size(inlier_loop_closures,1) * lc_outlier_prob / (1-lc_outlier_prob);
    
    while size(outlier_loop_closures,1) < num_outliers
        % randomly sample 2 vertices without replacement
        edge = randsample(n, 2);
        s = edge(1);
        t = edge(2);
        
        % make sure that the proposed edge does not already exist
        if ismember([s t], edges, 'rows')
            continue;
        end
        
        Rsrc = xtrue.R(:, (s-1)*d+1 : s*d);
        Rdst = xtrue.R(:, (t-1)*d+1 : t*d);
        tsrc = xtrue.t(:, s);
        tdst = xtrue.t(:, t);

        % get ground truth transformation
        Ts = [Rsrc        tsrc;
              zeros(1,d)     1];
        Tt = [Rdst        tdst;
              zeros(1,d)     1];
        dT = inv(Ts)*Tt;
        dR = dT(1:d,1:d);
        dt = dT(1:d,end);
        
        % Uniformly random relative translation
        error_t = -randtmax + 2*randtmax*rand(d,1);
        
        % Uniformly random relative rotation
        if d == 3
            error_R = randrot(3,1);
        else
            th = pi * rand;
            error_R = [cos(th), -sin(th); 
                       sin(th), cos(th)];
        end
        
        edge_idx = size(edges, 1) + 1;
        edges = [edges; s t];
        cell_t{edge_idx} = error_t;
        cell_R{edge_idx} = error_R;
        cell_tau{edge_idx} = tau;
        cell_kappa{edge_idx} = kappa;
        outlier_loop_closures = [outlier_loop_closures; s t];
    end
end

%% Return
new_measurements.edges = edges;
new_measurements.t = cell_t;
new_measurements.R = cell_R;
new_measurements.tau = cell_tau;
new_measurements.kappa = cell_kappa;

gt_info.inlier_loop_closures = inlier_loop_closures;
gt_info.outlier_loop_closures = outlier_loop_closures;

end