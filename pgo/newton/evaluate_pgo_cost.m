% function cost = evaluate_pgo_cost(measurements, R, t, options)
% Compute the cost of pose graph optimization.
% 
% Yulun Tian
function cost = evaluate_pgo_cost(measurements, R, t, options)
cost = 0;
d = size(measurements.t{1},1); 
assert(d == 3);
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
    ti = t(:,i);
    tj = t(:,j);
    kappa = measurements.kappa{k};
    tau = measurements.tau{k};
    Rij = measurements.R{k};
    tij = measurements.t{k};
    
    cost = cost + evaluate_relative_pose_cost(Ri, ti, Rj, tj, Rij, tij, kappa, tau, options);
end


end