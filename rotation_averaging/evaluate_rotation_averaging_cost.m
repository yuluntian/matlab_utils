% function cost = evaluate_rotation_averaging_cost(measurements, R, options, problem_data)
% Compute the cost of a rotation averaging problem
function cost = evaluate_rotation_averaging_cost(measurements, R, options, problem_data)
if nargin < 3
    options = struct;
end
if nargin < 4
    problem_data = struct;
end
if ~isfield(options, 'rotation_distance')
    options.rotation_distance = 'chordal';
end

% For chordal distance, if connection laplacian is also provided, use the
% faster method
if strcmp(options.rotation_distance, 'chordal') && isfield(problem_data, 'ConLap')
    cost = evaluate_rotation_averaging_cost_chordal(measurements, R, options, problem_data);
    return;
end

% Use default implementation. This might be slow if the problem has many
% edges.
cost = 0;
d = size(measurements.R{1},1); 
n = max(max(measurements.edges));
m = size(measurements.edges, 1);

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
    
    cost = cost + evaluate_relative_rotation_cost(Ri, Rj, Rij, kappa, options);
end


end

