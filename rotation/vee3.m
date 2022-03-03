% function w = vee3(W)
% Convert the element of so(3), i.e., skew-symmetric matrix to the
% corresponding 3 vector.
function w = vee3(W)

assert(size(W, 1) == 3);
assert(size(W, 2) == 3);
assert(issymmetric(W, 'skew'));

w = [W(3,2);
     W(1,3);
     W(2,1)];

end