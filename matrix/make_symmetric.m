% function [A, error] = make_symmetric(B)
% Create the symmetric version of the input matrix B
% Optionally compute the difference in Frobenieus norm
% between the original and symmetric version.
%
% Yulun Tian
function [A, diff] = make_symmetric(B)
A = 0.5 * (B + B');
if nargout > 1
    diff = norm(A - B, 'fro');
end
end