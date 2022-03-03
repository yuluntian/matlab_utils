clear; clc;

for trial = 1:1000
    w = randn(3,1);
    W = hat3(w);
    R1 = exp3(w);
    R2 = expm(W);
    assert(norm(R1-R2, 'fro') < 1e-6);
end

fprintf('Ok.\n');