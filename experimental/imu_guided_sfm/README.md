# IMU-Guided Essential Matrix Estimation

This experimental module provides a 5x-10x speedup in geometric verification by using IMU rotation data to reduce the Essential Matrix estimation from a 5-point problem to a 2-point problem.

## Implementation Details

- **Solver**: `solver_essential_two_point_imu.hpp` uses 2 point correspondences and a known relative rotation $R = R_2 R_1^T$ to solve for translation $t$.
- **Kernel**: `E_ACRobust_Imu.hpp` uses **Angular Error** (degrees/radians), making it completely independent of image resolution.

## Tuning Guide

| Threshold | Unit | Default | Tuning Advice |
| :--- | :--- | :--- | :--- |
| **Max Angular Error** | Degrees | `0.2` | **Resolution-independent.** Use `0.1` for high precision, `0.5` if IMU is noisy. |
| **RANSAC Iterations** | Count | `256` | Budget for finding the best translation direction. |
| **Min Inlier Count** | Count | `15` | Minimum consensus to accept a pair; higher prevents spurious matches. |

## Why it works
By using **Angular Error**, we measure the physical deviation of the light rays in 3D space. 

Traditional pixel-based thresholds (`4px`) are brittle because they change their physical meaning as resolution increases. A `0.2°` threshold remains constant whether the image is 4K or VGA, ensuring consistent filtering behavior across all sensor types and downsampling levels.

Additionally, by fixing the 3-DOF rotation using IMU, we only search a 2-DOF space (translation direction), which is significantly more robust against outliers and much faster to calculate than the standard 5-point algorithm.
