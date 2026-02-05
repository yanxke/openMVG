// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_ENGINE_HPP
#define OPENMVG_SFM_SFM_ENGINE_HPP

#include <string>

#include "openMVG/cameras/Camera_Common.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"

namespace openMVG {
namespace sfm {

/// Basic Reconstruction Engine.
/// Process Function handle the reconstruction.
class ReconstructionEngine
{
public:

  ReconstructionEngine
  (
    const SfM_Data & sfm_data,
    const std::string & soutDirectory
  )
  :sOut_directory_(soutDirectory),
    sfm_data_(sfm_data),
    intrinsic_refinement_options_(cameras::Intrinsic_Parameter_Type::ADJUST_ALL),
    extrinsic_refinement_options_(sfm::Extrinsic_Parameter_Type::ADJUST_ALL),
    b_use_motion_prior_(false),
    ba_rotation_steps_(0),
    ba_intrinsics_steps_(0),
    ba_final_steps_(0)
  {
  }

  virtual ~ReconstructionEngine() = default;

  virtual bool Process() = 0;

  cameras::Intrinsic_Parameter_Type Get_Intrinsics_Refinement_Type() const
  {
    return intrinsic_refinement_options_;
  }

  void Set_Intrinsics_Refinement_Type
  (
    cameras::Intrinsic_Parameter_Type rhs
  )
  {
    intrinsic_refinement_options_ = rhs;
  }

  Extrinsic_Parameter_Type Get_Extrinsics_Refinement_Type() const
  {
    return extrinsic_refinement_options_;
  }

  void Set_Extrinsics_Refinement_Type
  (
    sfm::Extrinsic_Parameter_Type rhs
  )
  {
    extrinsic_refinement_options_ = rhs;
  }

  void Set_Use_Motion_Prior
  (
    bool rhs
  )
  {
    b_use_motion_prior_ = rhs;
  }

  void Set_BA_Rotation_Steps(int steps)
  {
    ba_rotation_steps_ = steps;
  }

  int Get_BA_Rotation_Steps() const
  {
    return ba_rotation_steps_;
  }

  void Set_BA_Intrinsics_Steps(int steps)
  {
    ba_intrinsics_steps_ = steps;
  }

  int Get_BA_Intrinsics_Steps() const
  {
    return ba_intrinsics_steps_;
  }

  void Set_BA_Final_Steps(int steps)
  {
    ba_final_steps_ = steps;
  }

  int Get_BA_Final_Steps() const
  {
    return ba_final_steps_;
  }

  const SfM_Data & Get_SfM_Data() const {return sfm_data_;}

protected:
  std::string sOut_directory_; // Output path where outputs will be stored

  //-----
  //-- Reconstruction data
  //-----
  SfM_Data sfm_data_; // internal SfM_Data

  //-----
  //-- Reconstruction parameters
  //-----
  cameras::Intrinsic_Parameter_Type intrinsic_refinement_options_;
  sfm::Extrinsic_Parameter_Type extrinsic_refinement_options_;
  bool b_use_motion_prior_;

  // Bundle Adjustment step control (-1=disable, 0=use default, >0=limit steps)
  int ba_rotation_steps_;
  int ba_intrinsics_steps_;
  int ba_final_steps_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_ENGINE_HPP
