function [measurements, poses] = load_from_g2o(g2o_data_file)
% function [measurements, poses] = load_from_g2o(g2o_data_file)
% This function accepts as input a .g2o file containing the description of
% a 2D or 3D pose graph SLAM problem, and returns a MATLAB struct
% 'measurements' containing the description of the problem in the format 
% required by SE-Sync.  
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
% 'poses' consists of:
%
% vertices: (nx1)-dimensional vector encoding vertex id
% R: n-dimensional array of rotation matrices
% t: n-dimensional array of translation vectors
% 
% Author: Yulun Tian

fid = fopen(g2o_data_file, 'r');
edge_id = 0;
read_line = fgets(fid);  % Read the next line from the file

vertices = [];
pR = {};
pt = {};

while ischar(read_line)  % As long as this line is a valid character string
    
    token = strtok(read_line);
    
    if(strcmp(token, 'VERTEX_SE3:QUAT'))
        % 3D pose format:
        % VERTEX_SE3:QUAT id x y z qx qy qz qw
        C = textscan(read_line, '%s %d64  %f %f %f  %f %f %f %f');
        [name, id, x, y, z, qx, qy, qz, qw] = C{:};
        
        % Store the vertex id
        vertices = [vertices id + 1];  % MATLAB uses 1-based indexing
        
        % Store the translation vector
        pt{end+1} = [x, y, z]';
        
        q = [qw, qx, qy, qz]';
        q = q / norm(q);  % Make sure that this is properly normalized
        
        % Compute and store corresponding rotation matrix
        pR{end+1} = quat2rot(q);
        
    elseif(strcmp(token, 'VERTEX_SE2'))
        % 2D pose format:
        % VERTEX_SE2 ID x_meters y_meters yaw_radians
        C = textscan(read_line, '%s %d64  %f %f %f');
        [name, id, x, y, th] = C{:};
        
        % Store the vertex id
        vertices = [vertices id + 1];  % MATLAB uses 1-based indexing
        
        % Store the translation vector
        pt{end+1} = [x, y]';
       
        % Reconstruct and store the rotation
        pR{end+1} = [cos(th), -sin(th); 
                     sin(th), cos(th)];
    
    elseif(strcmp(token, 'EDGE_SE3:QUAT'))
        % 3D OBSERVATION
    
        edge_id = edge_id + 1;  % Increment the count for the number of edges
        
        
        % The g2o format specifies a 3D relative pose measurement in the
        % following form:
        
        % EDGE_SE3:QUAT id1 id2 dx dy dz dqx dqy dqz dqw
        % I11 I12 I13 I14 I15 I16
        %     I22 I23 I24 I25 I26
        %         I33 I34 I35 I36
        %             I44 I45 I46
        %                 I55 I56
        %                     I66
        C = textscan(read_line, '%s %d64 %d64 %f %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f   %f %f %f %f   %f %f %f    %f %f    %f');
        [name, id1, id2, dx, dy, dz, dqx, dqy, dqz, dqw, ...
            I11, I12, I13, I14, I15, I16, ...
            I22, I23, I24, I25, I26, ...
            I33, I34, I35, I36, ...
            I44, I45, I46, ...
            I55, I56, ...
            I66] = C{:};
        
        % Store the connectivity of this edge
        edges(edge_id, :) = [id1 + 1, id2 + 1];  % NB: .g2o uses 0-based indexing, whereas MATLAB uses 1-based indexing
        
        % Store the translational measurement
        t{edge_id} = [dx, dy, dz]';
        
        % Reconstruct quaternion for relative measurement
        q = [dqw, dqx, dqy, dqz]';
        q = q / norm(q);  % Make sure that this is properly normalized
        
        % Compute and store corresponding rotation matrix
        R{edge_id} = quat2rot(q);
        
        % Reconstruct the information matrix
        measurement_info = ...
            [I11, I12, I13, I14, I15, I16,
            I12, I22, I23, I24, I25, I26,
            I13, I23, I33, I34, I35, I36,
            I14, I24, I34, I44, I45, I46,
            I15, I25, I35, I45, I55, I56,
            I16, I26, I36, I46, I56, I66];
        
        % Compute and store the optimal (information-divergence-minimizing)
        % value of the parameter tau
        
        tau{edge_id} = 3 / trace(inv(measurement_info(1:3, 1:3)) );
        
        
        % Extract and store the optimal (information-divergence-minimizing)
        % value of the parameter kappa
        kappa{edge_id} = 3 / (2 *trace(inv(measurement_info(4:6, 4:6))));
    
    elseif(strcmp(token, 'EDGE_SE2'))
        % 2D OBSERVATION
        
        edge_id = edge_id + 1;
        
         % The g2o format specifies a 3D relative pose measurement in the
        % following form:
        
        % EDGE_SE2 id1 id2 dx dy dtheta, I11, I12, I13, I22, I23, I33
        C = textscan(read_line, '%s %d64 %d64 %f %f %f %f %f %f %f %f %f');
        [name, id1, id2, dx, dy, dth, I11, I12, I13, I22, I23, I33] = C{:};
        
        % Store the connectivity of this edge
        edges(edge_id, :) = [id1 + 1, id2 + 1];  % NB: .g2o uses 0-based indexing, whereas MATLAB uses 1-based indexing
        
        % Store the translational measurement
        t{edge_id} = [dx, dy]';
       
        % Reconstruct and store the rotational measurement
        R{edge_id} = [cos(dth), -sin(dth); 
                     sin(dth), cos(dth)];
        
        % Reconstruct the information matrix
        measurement_info = ...
            [I11, I12, I13;
            I12, I22, I23;
            I13, I23, I33];
        
        % Extract and store an outer approximation for the translational 
        % measurement precision
        tau{edge_id} = 2 / trace( inv(measurement_info(1:2, 1:2)));
        
        % Extract and store an outer approximation for the rotational
        % measurement precision
        kappa{edge_id} = I33;
       
    else
        error('Unknown name while reading g2o file: %s', token);
    end
    
    read_line = fgets(fid);
end

% Construct and return measurements struct
measurements.edges = edges;
measurements.R = R;
measurements.t = t;
measurements.kappa = kappa;
measurements.tau = tau;

% Construct and return pose struct
poses.vertices = vertices;
poses.R = pR;
poses.t = pt;

fclose(fid);
end

