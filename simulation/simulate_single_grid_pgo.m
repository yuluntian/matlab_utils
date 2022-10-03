% function [measurements, true_pose, gt_info] = simulate_single_grid_pgo(options)
% Yulun Tian
function [measurements, true_pose, gt_info] = simulate_single_grid_pgo(options)
if nargin < 1
    options = struct;
end
% standard deviation of translation measurements in meter
if ~isfield(options, 't_stddev')
    options.t_stddev = 0.1;
end
% standard deviation of rotation measurements in degree
if ~isfield(options, 'deg_stddev')
    options.deg_stddev = 3;
end
if ~isfield(options, 'num_rows')
    options.num_rows = 3;
end
if ~isfield(options, 'num_columns')
    options.num_columns = 3;
end
if ~isfield(options, 'num_floors')
    options.num_floors = 3;
end
if ~isfield(options, 'fov')
    options.fov = 1.1;
end
if ~isfield(options, 'lc_prob')
    options.lc_prob = 0.1;
end
if isfield(options, 'offset')
    offset = options.offset;
else
    offset = 0;
end
if ~isfield(options, 'random_rotation')
   random_rotation =  false;
else
   random_rotation =  options.random_rotation;
end
% SUPPORT FOR SIMULATING OUTLIER MEASUREMENTS
% Percentage of outliers within all loop closures
if isfield(options, 'lc_outlier_prob')
    lc_outlier_prob = options.lc_outlier_prob;
    assert(lc_outlier_prob >= 0 && lc_outlier_prob <= 1);
else
    lc_outlier_prob = 0;
end
if isfield(options, 'consecutive_odometry')
    consecutive_odometry = options.consecutive_odometry;
else
    consecutive_odometry = false;
end

num_rows = options.num_rows;
num_columns = options.num_columns;
num_floors = options.num_floors;
num_poses_per_floor = num_rows * num_columns;
num_poses = num_poses_per_floor * num_floors;
fov = options.fov;
lc_prob = options.lc_prob;
t_stddev = options.t_stddev; % standard deviation of translation measurements in meter
deg_stddev = options.deg_stddev; % standard deviation of rotation measurements in degree
rad_stddev = deg2rad(deg_stddev); % standard deviation of rotation measurements in radian
tau = 1/(t_stddev.^2); % translational measurement weights
kappa = (1/rad_stddev.^2)/2; % rotational measurement weights

% Support noiseless measurements
noiseless_t = false;
noiseless_R = false;
if t_stddev < 1e-14
    fprintf('Simulating noiseless translation for inlier measurements.\n');
    noiseless_t = true;
    tau = 100;  % use fixed weight in noiseless case
end
if deg_stddev < 1e-14
    fprintf('Simulating noiseless rotation for inlier measurements.\n');
    noiseless_R = true;
    kappa = 10000;  % use fixed weight in noiseless case
end

%% create ground truth poses
for floor = 1:num_floors
    for row = 1:num_rows
       for col = 1:num_columns
           
           % row major indexing
           if mod(row,2) == 1
               floor_index = (row-1)*num_columns + col;
           else
               floor_index = (row-1)*num_columns + (num_columns-col+1);
           end
           
           index = (floor-1)*num_poses_per_floor + floor_index;
           index_to_vertex{index} = [row; col; floor];
           vertex_to_index(row, col, floor) = index;
           
           % translation confined to be on the 3D grid
           true_pose.t{index} = offset + [row; col; floor];
           
       end
    end
end

% create ground truth rotation
for floor = 1:num_floors
    for row = 1:num_rows
       for col = 1:num_columns
           index = vertex_to_index(row, col, floor);
           
           % default euler angles
           yaw = 0;
           roll = 0;
           pitch = 0;
           
           if index < length(index_to_vertex)
               next_vertex = index_to_vertex{index+1};
               next_row = next_vertex(1);
               next_col = next_vertex(2);
               next_floor = next_vertex(3);
               if floor == next_floor
                   % if this pose is not the last one on this floor, set its
                   % yaw as its current heading
                   delta_row = next_row - row;
                   delta_col = next_col - col;
                   yaw = atan2(delta_col, delta_row);
               end
           end
           
           % simulate disturbances in euler angles
           if ~random_rotation
               yaw = yaw + normrnd(0, 0.1);
               roll = roll + normrnd(0, 0.1);
               pitch = pitch + normrnd(0, 0.1);
           else
               yaw = unifrnd(0,2*pi);
               roll = unifrnd(0,2*pi);
               pitch = unifrnd(0,2*pi);
           end
           
           % convert euler angles to rotation matrix 
           % rotation around Z (yaw), Y (pitch), X (roll)
           if index == 1
               true_pose.R{index} = eul2rotm([0 0 0], 'ZYX');
           else
               true_pose.R{index} = eul2rotm([yaw pitch roll], 'ZYX');
           end
       end
    end
end


Rtrue = [true_pose.R{:}];
ttrue = [true_pose.t{:}];

inlier_loop_closures = [];
outlier_loop_closures = [];

%% create odometry measurements within each floor

measurements.edges = [];

for floor = 1:num_floors
    index_offset = vertex_to_index(1,1,floor) - 1; 
    for e = 1:num_poses_per_floor-1
        s = index_offset + e; % source vertex id
        t = index_offset + e + 1; % destination vertex id
        
        ne = size(measurements.edges,1) + 1;
        measurements.edges(ne, 1) = s;
        measurements.edges(ne, 2) = t;
        
        % get ground truth transformation
        Ts = [true_pose.R{s} true_pose.t{s};
              zeros(1,3)     1];
        Tt = [true_pose.R{t} true_pose.t{t};
              zeros(1,3)     1];
        dT = inv(Ts)*Tt;
        dR = dT(1:3,1:3);
        dt = dT(1:3,end);
        
        % simulate noisy measurements
        error_t = normrnd(0, t_stddev, [3,1]);
        error_R = Langevin_sampler_3D(eye(3), kappa);
        if noiseless_t
            error_t = zeros(3,1);
        end
        if noiseless_R
            error_R = eye(3);
        end
        measurements.t{ne} = dt + error_t;
        measurements.R{ne} = dR * error_R;
        measurements.tau{ne} = tau;
        measurements.kappa{ne} = kappa;
        
    end
end

% create odometry between floors
for floor = 1:num_floors-1
    if ~consecutive_odometry
        s = vertex_to_index(1,1,floor);
        t = vertex_to_index(1,1,floor+1);
    else
        s = floor * num_rows * num_columns;  % last vertex in the current floor
        t = s + 1;
        vertex_src = index_to_vertex{s};
        vertex_dst = index_to_vertex{t};
        floor_src = vertex_src(3);
        floor_dst = vertex_dst(3);
        assert(floor_src == floor);
        assert(floor_dst == floor + 1);
    end
    
    ne = size(measurements.edges,1) + 1;
    measurements.edges(ne, 1) = s;
    measurements.edges(ne, 2) = t;

    % get ground truth transformation
    Ts = [true_pose.R{s} true_pose.t{s};
          zeros(1,3)     1];
    Tt = [true_pose.R{t} true_pose.t{t};
          zeros(1,3)     1];
    dT = inv(Ts)*Tt;
    dR = dT(1:3,1:3);
    dt = dT(1:3,end);

    % simulate noisy measurements
    error_t = normrnd(0, t_stddev, [3,1]);
    error_R = Langevin_sampler_3D(eye(3), kappa);
    if noiseless_t
        error_t = zeros(3,1);
    end
    if noiseless_R
        error_R = eye(3);
    end
    measurements.t{ne} = dt + error_t;
    measurements.R{ne} = dR * error_R;
    measurements.tau{ne} = tau;
    measurements.kappa{ne} = kappa;
    
end

slist = measurements.edges(:,1);
tlist = measurements.edges(:,2);
G = graph(slist, tlist);
Aodom = adjacency(G);


%% randomly create loop closure edges

for s = 1:num_poses
    for t = s+1:num_poses
        if (Aodom(s,t)==0)...
                && (norm(true_pose.t{s} - true_pose.t{t}) < fov)...
                && randsample([true false],1,true,[lc_prob,1-lc_prob])
            ne = size(measurements.edges,1) + 1;
            measurements.edges(ne,1) = s;
            measurements.edges(ne,2) = t;
            
            % get ground truth transformation
            Ts = [true_pose.R{s} true_pose.t{s};
                  zeros(1,3)     1];
            Tt = [true_pose.R{t} true_pose.t{t};
                  zeros(1,3)     1];
            dT = inv(Ts)*Tt;
            dR = dT(1:3,1:3);
            dt = dT(1:3,end);
            
            % determine if this measurement is an outlier
            is_outlier = randsample([true false],1,true,...
                                    [lc_outlier_prob, 1-lc_outlier_prob]);

            if ~is_outlier
                % simulate inlier noise
                error_t = normrnd(0, t_stddev, [3,1]);
                error_R = Langevin_sampler_3D(eye(3), kappa);
                if noiseless_t
                    error_t = zeros(3,1);
                end
                if noiseless_R
                    error_R = eye(3);
                end
                inlier_loop_closures = [inlier_loop_closures; s t];
            else
                % simulate outlier noise
                error_t = -10 + 20*rand(3,1);
                error_R = randrot(3,1);
                outlier_loop_closures = [outlier_loop_closures; s t];
            end
            
            if ~is_outlier
                measurements.t{ne} = dt + error_t;
                measurements.R{ne} = dR * error_R;
            else
                measurements.t{ne} = error_t;
                measurements.R{ne} = error_R;
            end
            
            % Measurement weights correspond to inlier models, 
            % Even for outlier measurements
            measurements.tau{ne} = tau;
            measurements.kappa{ne} = kappa;
        end
    end
end

gt_info = struct('inlier_loop_closures', inlier_loop_closures, ...
                 'outlier_loop_closures', outlier_loop_closures);

%% plot
% plot_poses(ttrue, Rtrue, measurements.edges);

end



