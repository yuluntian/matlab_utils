% Given two sets of positions, compute the RMSE absolute trajectory error
%
% Each set of translation is a d-by-n matrix, where n is the number of
% points.
% 
% This function assumes that the input translations are already aligned in
% a common frame.
%
% Definition of ATE:
% Sturm et al. "A Benchmark for the Evaluation of RGB-D SLAM Systems"
% 
% Yulun Tian 

function ATE = compute_ATE(translations1, translations2)

d = size(translations1, 1);
n = size(translations1, 2);
assert(size(translations2, 1) == d);
assert(size(translations2, 2) == n);

% Compute the 2-norm for each column
distances = vecnorm(translations1 - translations2, 2, 1);
distances_squared = distances.*distances;

% Compute average
ATE = sqrt(mean(distances_squared));

end