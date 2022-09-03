clear; clc;

TESTNAME = 'test_SOR';

fprintf('Begin %s...\n', TESTNAME);

%% size of the linear system
n = 500;
num_blocks = 5;
poses_per_block = floor(n / num_blocks);
for b = 1:num_blocks
   if b < num_blocks
        blocks{b} = (b-1)*poses_per_block+1 : b*poses_per_block;
   else
        blocks{b} = (b-1)*poses_per_block+1 : n;
   end
end

%% Generate random positive definite linear system Ax = b
A = randn(n, n);
A = A * A' + 0.1 * eye(n); % make sure A is positive definite
b = randn(n, 1);
x0 = zeros(n, 1);
SOR_opts = struct;
SOR_opts.blocks = blocks;
SOR_opts.relchange_tol = 1e-8;
SOR_opts.w = 1.9;
SOR_opts.maxiter = 10000000;
SOR_opts.verbose = true;
SOR_opts.check_spd = true;
SOR_opts.check_spectral_radius = true;

%% Test!
xopt = A \ b;
[x, SOR_info] = SOR(A, b, x0, SOR_opts);
diff = norm(x - xopt);
assert(diff < 1e-4, ...
       'Result not consistent with centralized solution! Diff = %g', diff);
   
%% Exit
fprintf('%s OK.\n', TESTNAME);