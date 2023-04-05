% function t_dst = transform_translations(t_src, T_dst_src)
% Transform a set of translations from the input 'src' frame to the output 
% 'dst' frame by applying the known transformation.
%
% Each set of translation is a d-by-n matrix t = [t1 ... tn]
% 
% Yulun Tian 
function t_dst = transform_translations(t_src, T_dst_src)

d = size(t_src, 1);
n = size(t_src, 2);
assert(size(T_dst_src.t, 1) == d);
if isfield(T_dst_src, 's')
    warning('Applying similarity transform! (scale=%.2f)', T_dst_src.s);
    s = T_dst_src.s;
else
    s = 1.0;
end

t_dst = s * T_dst_src.R * t_src + repmat(T_dst_src.t, 1, n);

end