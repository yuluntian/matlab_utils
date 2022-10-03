% function [measurements, true_pose, blocks] = simulate_multi_grid_pgo(sim_opts)
% Yulun Tian
function [measurements, true_pose, blocks] = simulate_multi_grid_pgo(sim_opts)
num_robot_rows = sim_opts.num_robot_rows;
num_robot_columns = sim_opts.num_robot_columns;
num_rows = sim_opts.num_rows; % height of a single grid
num_columns = sim_opts.num_columns; % width of a single grid
num_floors  = sim_opts.num_floors; % height of a single grid
inter_lc_prob = sim_opts.inter_lc_prob; % probability of inter-robot loop closure
fov = sim_opts.fov; % field of view, within which loop closure is sampled
t_stddev = sim_opts.t_stddev; % translational measurement standard deviation (meter)
deg_stddev = sim_opts.deg_stddev; % rotational measurement standard deviation (degree)


rad_stddev = deg2rad(deg_stddev); % standard deviation of rotation measurements in radian
tau = 1/(t_stddev.^2); % translational measurement weights
kappa = (1/rad_stddev.^2)/2; % rotational measurement weights

%% simulate grid of each robot
num_robots = num_robot_rows * num_robot_columns;
num_pose_per_robot = num_rows * num_columns * num_floors;
meas_struct = cell(1, num_robots);
tp_struct = cell(1, num_robots);
blocks = cell(1, num_robots);

for robot_col = 1:num_robot_columns
    for robot_row = 1:num_robot_rows
       robot_index = (robot_row-1)*num_robot_columns + robot_col;
       sim_opts.offset = [(robot_col-1)*num_columns; (robot_row-1)*num_rows; 0];
       [meas, tp, ~] = simulate_single_grid_pgo(sim_opts);
       meas_struct{robot_index} = meas;
       tp_struct{robot_index} = tp;
       blocks{robot_index} = (robot_index-1)*num_pose_per_robot+1 : robot_index*num_pose_per_robot;
    end
end

%% merge true pose and intra-robot edges
for r = 1:num_robots
   for i = 1:num_pose_per_robot
      tp = tp_struct{r};
      merged_id = (r-1)*num_pose_per_robot + i;
      true_pose.t{merged_id} = tp.t{i};
      true_pose.R{merged_id} = tp.R{i};
   end
end

measurements.edges = [];
for r = 1:num_robots
    meas = meas_struct{r};
    for e = 1:size(meas.edges,1)
        edge_id = size(measurements.edges,1)+1;
        s = meas.edges(e,1) + (r-1)*num_pose_per_robot;
        t = meas.edges(e,2) + (r-1)*num_pose_per_robot;
        measurements.edges(edge_id,1) = s;
        measurements.edges(edge_id,2) = t;
        measurements.R{edge_id} = meas.R{e};
        measurements.t{edge_id} = meas.t{e};
        measurements.kappa{edge_id} = meas.kappa{e};
        measurements.tau{edge_id} = meas.tau{e};
    end
end

%% create inter robot loop closures
for r1 = 1:num_robots
   for r2 = r1+1:num_robots
      inter_count = 0;
      for s = blocks{r1}
         for t = blocks{r2}
            if norm(true_pose.t{s} - true_pose.t{t}) < fov
                % ensure that there is at least one inter robot loop
                % closure so that graph is connected
                if inter_count == 0 || randsample([true false],1,true,[inter_lc_prob,1-inter_lc_prob])
                    inter_count = inter_count + 1;
                    edge_id = size(measurements.edges,1) + 1;
                    measurements.edges(edge_id,1) = s;
                    measurements.edges(edge_id,2) = t;
                    
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
                    measurements.t{edge_id} = dt + error_t;
                    measurements.R{edge_id} = dR * error_R;
                    measurements.tau{edge_id} = tau;
                    measurements.kappa{edge_id} = kappa;
                end
            end
         end
      end
   end
end


end
