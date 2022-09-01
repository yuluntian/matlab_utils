clear; clc; close all;

% Test same block sizes
p = 2;
num_blocks = 4;
n = p * num_blocks;
A1 = randn(n,n);
D1 = extract_blkdiag_matrix(A1, p);

% Test different block sizes
bdims = [1 5 3];
n = sum(bdims);
A2 = randn(n,n);
D2 = extract_blkdiag_matrix(A2, bdims);