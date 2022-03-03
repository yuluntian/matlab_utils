clear; clc;
for trial = 1:100
    w = randn(3,1);
    W = hat3(w);
    w2 = vee3(W);
    assert(norm(w2-w) < 1e-7);
end
fprintf('Ok.\n');