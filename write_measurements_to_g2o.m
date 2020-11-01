% Write a measurement struct to standard g2o format
% Input measurement struct format:
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
% Yulun Tian
function write_measurements_to_g2o(measurements, filename)
fout = fopen(filename, 'w'); % discard existing content if any

d = size(measurements.t{1},1);
assert(d == 3, 'Only support 3D measurements!');

% iterate over edges
for k = 1:size(measurements.edges,1)
    vertex_src = measurements.edges(k,1) - 1; % g2o uses 0-based indexing
    vertex_dst = measurements.edges(k,2) - 1;
    R = measurements.R{k};
    t = measurements.t{k};
    kappa = measurements.kappa{k};
    tau = measurements.tau{k};
    
    quat = rotm2quat(R);
    qw = quat(1); qx = quat(2); qy = quat(3); qz = quat(4);
    tx = t(1); ty = t(2); tz = t(3);
    
    I = zeros(6);
    I(1:3,1:3) = tau * eye(3);
    I(4:6,4:6) = 2 * kappa * eye(3);
    
    % The g2o format specifies a 3D relative pose measurement in the
        % following form:
        
    % EDGE_SE3:QUAT id1 id2 dx dy dz dqx dqy dqz dqw
    % I11 I12 I13 I14 I15 I16
    %     I22 I23 I24 I25 I26
    %         I33 I34 I35 I36
    %             I44 I45 I46
    %                 I55 I56
    %                     I66
    
    g2oLineSpec = 'EDGE_SE3:QUAT %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n';
    fprintf(fout, g2oLineSpec, ...
            vertex_src, vertex_dst, ...
            tx, ty, tz, ...
            qx, qy, qz, qw, ...
            I(1,1), I(1,2), I(1,3), I(1,4), I(1,5), I(1,6), ...
                    I(2,2), I(2,3), I(2,4), I(2,5), I(2,6), ...
                            I(3,3), I(3,4), I(3,5), I(3,6), ...
                                    I(4,4), I(4,5), I(4,6), ...
                                            I(5,5), I(5,6), ...
                                                    I(6,6));
end

fprintf('Wrote %i edges to %s \n', size(measurements.edges,1), filename);
end
