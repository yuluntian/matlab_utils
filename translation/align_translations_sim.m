% function [T_dst_src, translations_src_aligned, info] = align_translations_sim(translations_src, translations_dst)
% Given two sets of points in different global frames, find the similarity
% transformation (s, R, t) that transforms points in the src frame to points in the
% dst frame and minimizes the squared error.
% \sum || s R p_src + t - p_dst||^2
% Reference: http://graphics.stanford.edu/~smr/ICP/comparison/horn-hilden-orientation-josa88.pdf
%
% Each set of translation is a d-by-n matrix, where n is the number of
% points.
% 
% Yulun Tian 

function [T_dst_src, translations_src_aligned, info] = align_translations_sim(translations_src, translations_dst)

d = size(translations_src, 1);
n = size(translations_src, 2);
assert(d == 2 || d == 3);
assert(size(translations_dst, 1) == d);
assert(size(translations_dst, 2) == n);

cent_src = mean(translations_src, 2);
cent_dst = mean(translations_dst, 2);

t_src_centered = translations_src - repmat(cent_src, 1, n);
t_dst_centered = translations_dst - repmat(cent_dst, 1, n);

% recover scale
scale = sqrt(norm(t_dst_centered, 'fro').^2 / norm(t_src_centered, 'fro').^2);

% recover rotation
M = t_dst_centered * t_src_centered';
R = project_to_rotation_matrix(M);

% recover translation
t = cent_dst - scale * R * cent_src;

T_dst_src = struct;
T_dst_src.s = scale;
T_dst_src.R = R;
T_dst_src.t = t;

% return aligned translations
translations_src_aligned = scale * R * translations_src + repmat(t, 1, n);

% no additional info for now
info = struct;

end

