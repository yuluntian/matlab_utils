clear; clc;

for trial = 1:10000
    % Sample two rotations
    R1 = randrot(3);
    R2 = randrot(3);
    d = norm(R1-R2, 'fro');
    th = chordal_to_angular(d);
    d2 = angular_to_chordal(th);
    assert(abs(d-d2) < 1e-8);
end

% 2D 
for trial = 1:10000
    th1 = 2*pi*rand - pi;
    th2 = 2*pi*rand - pi;
    R1 = euler_to_rotm2d(th1);
    R2 = euler_to_rotm2d(th2);
    dchr = norm(R1 - R2, 'fro');
    dang = norm(log2d(R1' * R2));
    dchr2 = angular_to_chordal(dang);
    assert(abs(dchr - dchr2) < 1e-10);
end

fprintf('Ok.\n');