clear; clc; close all;

bdims = [1 4 3];
n = sum(bdims);
A = randn(n,n);
D = extract_blkdiag_matrix(A, bdims);