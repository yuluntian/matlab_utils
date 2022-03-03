% function R = exp3(w)
% Compute the exponential map for 3d rotation
% Input is either a 3-vector or the corresponding 3-by-3 skew symmetric
% matrix.
% Reference: eq (7.23) of Barfoot textbook
% 
% Yulun Tian
function R = exp3(x)
assert(size(x,1) == 3);

if size(x,2) == 1
    w = x;
elseif size(x,2) == 3
    w = vee3(x);
else
    error('Input is neither a 3 vector of a 3-by-3 skew symmetric matrix.');
end

th = norm(w);
a = w / th;

R = cos(th)*eye(3) + (1-cos(th)) * (a * a') + sin(th) * hat3(a);

end