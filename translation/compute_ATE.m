% function ATE = compute_ATE(translations1, translations2, options)
% Given two sets of positions, compute the RMSE absolute trajectory error.
%
% Each set of translation is a d-by-n matrix, where n is the number of
% points.
%
% options.align: if set to true, will first align the input positions in
% the global frame before computing ATE.
%
% Definition of ATE:
% Sturm et al. "A Benchmark for the Evaluation of RGB-D SLAM Systems"
% 
% Yulun Tian 

function ATE = compute_ATE(translations1, translations2, options)
if nargin < 3
    options = struct;
end
if ~isfield(options, 'align')
    options.align = true;
end

d = size(translations1, 1);
n = size(translations1, 2);
assert(size(translations2, 1) == d);
assert(size(translations2, 2) == n);

% Align translations
if options.align
    [~, translations1] = align_translations(translations1, translations2);
end

% Compute the 2-norm for each column
distances = vecnorm(translations1 - translations2, 2, 1);
distances_squared = distances.*distances;

% Compute average
ATE = sqrt(mean(distances_squared));

end