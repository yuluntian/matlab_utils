# MATLAB Utilities for 3D geometry and optimization 

A MATLAB library for pose graph optimization (PGO), rotation averaging, and related geometric computations on manifolds.

## Overview

This library provides tools for solving optimization problems in robotics and computer vision, particularly focusing on:
- **Pose Graph Optimization (PGO)**: Optimizing robot trajectories from relative pose measurements
- **Rotation Averaging**: Estimating absolute rotations from relative rotation measurements
- **Geometric Computations**: Lie group operations, manifold optimization, and matrix utilities

## Setup

Run `setup.m` to initialize the library:
```matlab
setup
```

This will:
- Import Manopt (if available) for manifold optimization
- Add all subdirectories to the MATLAB path

**Dependencies:**
- MATLAB (tested with recent versions)
- [Manopt](https://www.manopt.org/) (optional but recommended for optimization functions)

## Directory Structure

### Core Modules

#### `pgo/` - Pose Graph Optimization
Solve full pose graph optimization problems.

#### `rotation_averaging/` - Rotation Averaging
Estimate absolute rotations from relative rotation measurements.

#### `rotation/` - Rotation Utilities
Low-level rotation operations and conversions.

#### `translation/` - Translation Utilities
Operations for translation estimation and alignment.

#### `pose/` - Pose Utilities
Combined rotation and translation operations.

#### `simulation/` - Problem Generation
Generate synthetic datasets for testing and benchmarking.

### Supporting Modules

#### `matrix/` - Matrix Operations
Specialized matrix operations for graph and optimization problems.

#### `linear_solver/` - Linear System Solvers
Efficient solvers for structured linear systems.

#### `g2o/` - G2O File I/O
Interface with g2o format.

#### `gtsam/` - GTSAM Utilities
Helper functions for GTSAM integration.

#### `misc/` - Miscellaneous
General utility functions.

## Usage Examples

### Example 1: Simulate and Solve a PGO Problem

```matlab
% Generate synthetic PGO problem
options.num_rows = 5;
options.num_columns = 5;
options.deg_stddev = 3;  % 3 degrees rotation noise
options.t_stddev = 0.1;  % 0.1m translation noise
[measurements, true_pose] = simulate_single_grid_pgo(options);

% Solve using Gauss-Newton
R_init = eye(3, 3*numel(measurements.vertices));
t_init = zeros(3, numel(measurements.vertices));
[R_opt, t_opt] = pgo_gauss_newton(measurements, R_init, t_init);

% Evaluate cost
cost = evaluate_pgo_cost(measurements, R_opt, t_opt);
fprintf('Final cost: %.4f\n', cost);
```

### Example 2: Rotation Averaging

```matlab
% Extract rotation measurements
rot_measurements.edges = measurements.edges;
rot_measurements.R = measurements.R;
rot_measurements.kappa = measurements.kappa;

% Solve rotation averaging
[R_avg, info] = rotation_averaging_gauss_newton(rot_measurements);

% Compute error
cost = evaluate_rotation_averaging_cost(rot_measurements, R_avg);
```

### Example 3: Load from G2O File

```matlab
% Load dataset
[measurements, poses_gt] = load_from_g2o('dataset.g2o');

% Solve
[R_opt, t_opt] = pgo_gauss_newton(measurements, poses_gt.R, poses_gt.t);

% Save results
write_to_g2o('output.g2o', measurements, R_opt, t_opt);
```
## Author

Yulun Tian (yulunt@umich.edu)

## License

See individual files for licensing information.

## Acknowledgements
This repo built on the awesome implementations from the following prior work:

- Manopt: [www.manopt.org](https://www.manopt.org/)
- g2o format: [github.com/RainerKuemmerle/g2o](https://github.com/RainerKuemmerle/g2o)
- SE-Sync: Rosen et al., "SE-Sync: A certifiably correct algorithm for synchronization over the special Euclidean group" [github.com/david-m-rosen/SE-Sync](https://github.com/david-m-rosen/SE-Sync)
