% function [T_dst_src, translations_src_aligned] = align_translations(translations_src, translations_dst)
% Given two sets of points in different global frames, find the rigid-body
% transformation that transform points in the src frame to points in the
% dst frame.
%
% Each set of translation is a d-by-n matrix, where n is the number of
% points.
% 
% Example reference: 
% https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
% 
% Yulun Tian 

function [T_dst_src, translations_src_aligned] = align_translations(translations_src, translations_dst)

d = size(translations_src, 1);
n = size(translations_src, 2);
assert(d == 2 || d == 3);
assert(size(translations_dst, 1) == d);
assert(size(translations_dst, 2) == n);

cent_src = mean(translations_src, 2);
cent_dst = mean(translations_dst, 2);

H = 0;
for i = 1:n
    v_src = translations_src(:,i) - cent_src;
    v_dst = translations_dst(:,i) - cent_dst;
    H = H + v_src*v_dst';
end

[U,~,V] = svd(H);
reflector = eye(d);
reflector(end,end) = det(V*U');

% recover best aligning rotation and translation
R_dst_src = V * reflector * U';
t_dst_src = cent_dst - R_dst_src*cent_src;
T_dst_src.R = R_dst_src;
T_dst_src.t = t_dst_src;

% return aligned translations
translations_src_aligned = R_dst_src * translations_src + repmat(t_dst_src, 1, n);
end

