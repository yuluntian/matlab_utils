% function print_if(condition, string, varargin)
% Only print the input string to command line if the input condition
% evaluates to true.
%
% Yulun Tian.
function print_if(condition, string, varargin)
if condition
    fprintf(string, varargin{:});
end
end