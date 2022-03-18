clear; clc;

for trial = 1:1000
    th1 = pi * rand;
    R = euler_to_rotm2d(th1);
    th2 = rotm2d_to_euler(R);
    assert(abs(th1 - th2) < 1e-8);
end

fprintf('Ok.\n');