% function w = vee2(W)
% Convert the element of so(2), i.e., skew-symmetric matrix to the
% corresponding scalar.
function w = vee2(W)

assert(size(W, 1) == 2);
assert(size(W, 2) == 2);
%assert(issymmetric(W, 'skew'));
ss_err = norm(W+W', 'fro');
if ss_err > 1e-6
    warning('Input is not skew-symmetric!');
end

w = W(2,1);

end