clear; clc;

for trial = 1:10000
    % Sample two rotations
    R1 = randrot(3);
    R2 = randrot(3);
    d = norm(R1-R2, 'fro');
    th = chordal2angular(d);
    d2 = angular2chordal(th);
    assert(abs(d-d2) < 1e-8);
end

fprintf('Ok.\n');