/*
 * Integration Snippet for main_GeometricFilterWithPriors.cpp
 * 
 * 1. Add ESSENTIAL_MATRIX_IMU to EGeometricModel
 * 2. Add IMU parsing helpers (logic from sfm_global_engine_relative_motions.cpp)
 * 3. Add to switch statement:
 */

case ESSENTIAL_MATRIX_IMU:
{
  filter_ptr->Robust_model_estimation(
      // 0.2 is the angular precision in DEGREES (resolution-independent)
      GeometricFilter_EMatrix_AC_Imu( 0.2, 256, &map_imu_rotations, map_PutativeMatches.size() ),
      map_PutativeMatches,
      bGuided_matching,
      d_distance_ratio,
      &progress );
  map_GeometricMatches = filter_ptr->Get_geometric_matches();
}
break;
