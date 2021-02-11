% Given a set of translations, plot as a trajectory in the input axes
%
% The set of translation is a d-by-n matrix, where n is the number of
% points.
% 
% This function can optinally take in options.R and options.t to apply a
% global transformation of the points before plotting.
% 
% Yulun Tian 
function line = plot_translations(ax, translations, options)

d = size(translations, 1);
n = size(translations, 2);
assert(d == 2 || d == 3);

if nargin < 3
    options = struct;
end

% Check if we need to apply a global transformation before plotting
if(isfield(options, 'R'))
   assert(isfield(options, 't'),  'options.t must be privided!');
   translations = options.R * translations + repmat(options.t, 1, n);
end

if d == 2
   line = plot(ax, ...
               translations(1,:), ...
               translations(2,:));
else
   line = plot3(ax, ...
                translations(1,:), ...
                translations(2,:), ...
                translations(3,:));
end

% Pass in optional setups
if isfield(options, 'Color')
    line.Color = options.Color;
end

if isfield(options, 'LineWidth')
    line.LineWidth = options.LineWidth;
end

if isfield(options, 'LineStyle')
    line.LineStyle = options.LineStyle;
end

if isfield(options, 'DisplayName')
    line.DisplayName = options.DisplayName;
end
end