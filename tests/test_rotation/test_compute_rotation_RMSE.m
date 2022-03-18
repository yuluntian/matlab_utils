clear; clc; 

for trial = 1:100
    R1 = reshape(randrot(3,10), 3, []);
    R2 = R1;
    error_deg = compute_rotation_RMSE(R1,R2);
    assert(error_deg < 1e-6);
end

fprintf('Ok.\n');