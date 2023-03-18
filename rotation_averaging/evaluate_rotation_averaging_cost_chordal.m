% function cost = evaluate_rotation_averaging_cost_chordal(measurements, R, options, problem_data)
% Compute the rotation averaging cost function under the squared chordal
% distance.
%
% Yulun Tian
function cost = evaluate_rotation_averaging_cost_chordal(measurements, R, options, problem_data)
    assert(strcmp(options.rotation_distance, 'chordal'));
    assert(isfield(problem_data, 'ConLap'));
    cost = 0.5 * trace(R * (problem_data.ConLap * R'));
end