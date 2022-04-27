% function cost = evaluate_rotation_averaging_cost(measurements, R, options)
% Compute the cost of a rotation averaging problem
function cost = evaluate_rotation_averaging_cost(measurements, R, options)
if nargin < 3
    options = struct;
end
if ~isfield(options, 'rotation_distance')
    options.rotation_distance = 'chordal';
end

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