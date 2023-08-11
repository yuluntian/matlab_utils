% function T = combine_rotations_and_translations(R, t)
% Given R = [R1, ... Rn] and t = [t1, ..., tn]
% Return T = [R1 t1 ... Rn tn]
%
% Yulun Tian
function T = combine_rotations_and_translations(R, t)

d = size(t, 1);
n = size(t, 2);
assert(size(R,1) == d)
assert(size(R,2) == d*n)

T = [];
for i = 1:n
    Ri = R(:, (i-1)*d+1 : i*d);
    ti = t(:, i);
    T = [T Ri ti];
end

end