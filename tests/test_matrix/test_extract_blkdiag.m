clear; clc; close all;

p = 2;
num_blocks = 4;
n = p * num_blocks;
A = randn(n,n);
D = extract_blkdiag_matrix(A, p);