clear; clc;

w = 0;
W = zeros(2,2);
R1 = exp2(w);
R2 = expm(W);
assert(norm(R1-R2, 'fro') < 1e-6);

for trial = 1:1000
    w1 = randn(1,1);
    w2 = log2d(exp2(w1));
    assert(abs(w1-w2) < 1e-6);
end

fprintf('Ok.\n');