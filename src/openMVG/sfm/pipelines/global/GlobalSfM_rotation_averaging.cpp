// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/pipelines/global/GlobalSfM_rotation_averaging.hpp"

#include "openMVG/graph/graph.hpp"
#include "openMVG/multiview/rotation_averaging.hpp"
#include "openMVG/sfm/sfm_filters.hpp"
#include "openMVG/sfm/pipelines/global/sfm_global_reindex.hpp"
#include "openMVG/stl/stlMap.hpp"
#include "openMVG/system/logger.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/tracks/union_find.hpp"

#include "third_party/histogram/histogram.hpp"

#include <algorithm>
#include <map>
#include <set>

namespace openMVG{
namespace sfm{

using namespace openMVG::rotation_averaging;

Pair_Set GlobalSfM_Rotation_AveragingSolver::GetUsedPairs() const
{
  return used_pairs;
}

bool GlobalSfM_Rotation_AveragingSolver::Run(
  ERotationAveragingMethod eRotationAveragingMethod,
  ERelativeRotationInferenceMethod eRelativeRotationInferenceMethod,
  const RelativeRotations & relativeRot_In,
  Hash_Map<IndexT, Mat3> & map_globalR,
  IndexT fixed_pose_id,
  const Hash_Map<IndexT, double> * pose_timestamps,
  const Hash_Map<IndexT, std::string> * pose_img_names
) const
{
  RelativeRotations relativeRotations = relativeRot_In;
  // We work on a copy, since inference can remove some relative motions

  switch (eRelativeRotationInferenceMethod)
  {
    case TRIPLET_ROTATION_INFERENCE_NONE:
    break;
    case TRIPLET_ROTATION_INFERENCE_COMPOSITION_ERROR:
    {
      //-------------------
      // Triplet inference (test over the composition error)
      //-------------------
      Pair_Set pairs = getPairs(relativeRotations);
      std::vector<graph::Triplet> vec_triplets = graph::TripletListing(pairs);

      //-- Rejection triplet that are 'not' identity rotation (error to identity > 5°)
      TripletRotationRejection(5.0f, vec_triplets, relativeRotations, pose_timestamps, pose_img_names);

      pairs = getPairs(relativeRotations);
      const std::set<IndexT> set_remainingIds = graph::CleanGraph_KeepLargestBiEdge_Nodes<Pair_Set, IndexT>(pairs);
      if (set_remainingIds.empty())
        return false;
      KeepOnlyReferencedElement(set_remainingIds, relativeRotations);
    }
  break;
    default:
    OPENMVG_LOG_ERROR
      << "Unknown relative rotation inference method: "
      << (int) eRelativeRotationInferenceMethod;
  }

  // Compute contiguous index (mapping between sparse index and contiguous index)
  //  from ranging in [min(Id), max(Id)] to  [0, nbCam]

  const Pair_Set pairs = getPairs(relativeRotations);
  Hash_Map<IndexT, IndexT> reindexForward, reindexBackward;
  reindex(pairs, reindexForward, reindexBackward);

  for (RelativeRotations::iterator iter = relativeRotations.begin();  iter != relativeRotations.end(); ++iter)
  {
    RelativeRotation & rel = *iter;
    rel.i = reindexForward[rel.i];
    rel.j = reindexForward[rel.j];
  }

  if (fixed_pose_id != UndefinedIndexT)
  {
    const auto it = reindexForward.find(fixed_pose_id);
    if (it != reindexForward.end() && it->second != 0)
    {
      const IndexT fixed_idx = it->second;
      const IndexT old0 = reindexBackward[0];

      reindexForward[fixed_pose_id] = 0;
      reindexForward[old0] = fixed_idx;
      reindexBackward[0] = fixed_pose_id;
      reindexBackward[fixed_idx] = old0;

      for (RelativeRotations::iterator iter = relativeRotations.begin();  iter != relativeRotations.end(); ++iter)
      {
        RelativeRotation & rel = *iter;
        if (rel.i == 0) rel.i = fixed_idx;
        else if (rel.i == fixed_idx) rel.i = 0;
        if (rel.j == 0) rel.j = fixed_idx;
        else if (rel.j == fixed_idx) rel.j = 0;
      }
    }
  }

  //- B. solve global rotation computation
  bool bSuccess = false;
  std::vector<Mat3> vec_globalR(reindexForward.size());
  switch (eRotationAveragingMethod)
  {
    case ROTATION_AVERAGING_L2:
    {
      //- Solve the global rotation estimation problem:
      bSuccess = rotation_averaging::l2::L2RotationAveraging(
        reindexForward.size(),
        relativeRotations,
        vec_globalR);
      //- Non linear refinement of the global rotations
      if (bSuccess)
        bSuccess = rotation_averaging::l2::L2RotationAveraging_Refine(
          relativeRotations,
          vec_globalR);

      // save kept pairs (restore original pose indices using the backward reindexing)
      for (RelativeRotations::iterator iter = relativeRotations.begin();  iter != relativeRotations.end(); ++iter)
      {
        RelativeRotation & rel = *iter;
        rel.i = reindexBackward[rel.i];
        rel.j = reindexBackward[rel.j];
      }
      used_pairs = getPairs(relativeRotations);
    }
    break;
    case ROTATION_AVERAGING_L1:
    {
      using namespace openMVG::rotation_averaging::l1;

      //- Solve the global rotation estimation problem:
      const size_t nMainViewID = 0; //arbitrary choice
      std::vector<bool> vec_inliers;
      bSuccess = rotation_averaging::l1::GlobalRotationsRobust(
        relativeRotations, vec_globalR, nMainViewID, 0.0f, &vec_inliers);

      // Suppress verbose inlier mask logging.

      // save kept pairs (restore original pose indices using the backward reindexing)
      for (size_t i = 0; i < vec_inliers.size(); ++i)
      {
        if (vec_inliers[i])
        {
          used_pairs.insert(
            {reindexBackward[relativeRotations[i].i],
            reindexBackward[relativeRotations[i].j]});
        }
      }
    }
    break;
    default:
    OPENMVG_LOG_ERROR << "Unknown rotation averaging method: " << (int) eRotationAveragingMethod;
  }

  if (bSuccess)
  {
    //-- Setup the averaged rotations
    for (size_t i = 0; i < vec_globalR.size(); ++i)  {
      map_globalR[reindexBackward[i]] = vec_globalR[i];
    }
  }
  else {
    OPENMVG_LOG_ERROR << "Global rotation solving failed.";
  }

  return bSuccess;
}

void GlobalSfM_Rotation_AveragingSolver::TripletRotationRejection(
  const double max_angular_error,
  std::vector<graph::Triplet> & vec_triplets,
  RelativeRotations & relativeRotations,
  const Hash_Map<IndexT, double> * pose_timestamps,
  const Hash_Map<IndexT, std::string> * pose_img_names) const
{
  const size_t triplets_before_count = vec_triplets.size();
  const size_t edges_start_count = relativeRotations.size();

  RelativeRotations_map map_relatives = getMap(relativeRotations);
  RelativeRotations_map map_relatives_validated;

  struct TripletWithErr {
    graph::Triplet triplet;
    float err;
  };
  std::vector<TripletWithErr> all_triplets;
  all_triplets.reserve(vec_triplets.size());

  //--
  // 1. ROTATION ERROR COMPUTATION
  //--

  std::vector<float> vec_errToIdentityPerTriplet;
  vec_errToIdentityPerTriplet.reserve(vec_triplets.size());
  
  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graph::Triplet & triplet = vec_triplets[i];
    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;

    const Pair ij(I,J), ji(J,I);
    const Mat3 RIJ = (map_relatives.count(ij)) ?
      map_relatives.at(ij).Rij : Mat3(map_relatives.at(ji).Rij.transpose());

    const Pair jk(J,K), kj(K,J);
    const Mat3 RJK = (map_relatives.count(jk)) ?
      map_relatives.at(jk).Rij : Mat3(map_relatives.at(kj).Rij.transpose());

    const Pair ki(K,I), ik(I,K);
    const Mat3 RKI = (map_relatives.count(ki)) ?
      map_relatives.at(ki).Rij : Mat3(map_relatives.at(ik).Rij.transpose());

    const Mat3 Rot_To_Identity = RIJ * RJK * RKI;
    const float angularErrorDegree = static_cast<float>(R2D(getRotationMagnitude(Rot_To_Identity)));
    vec_errToIdentityPerTriplet.push_back(angularErrorDegree);

    if (angularErrorDegree <= 15.0f) {
      all_triplets.push_back({triplet, angularErrorDegree});
    }
  }

  //--
  // 2. CONNECTIVITY RESCUE PASS
  //--

  std::vector<graph::Triplet> vec_triplets_validated;
  std::set<size_t> validated_triplet_indices;

  if (pose_timestamps != nullptr && !all_triplets.empty())
  {
    std::set<IndexT> all_poses;
    for (const auto & tr : all_triplets) {
      all_poses.insert(tr.triplet.i);
      all_poses.insert(tr.triplet.j);
      all_poses.insert(tr.triplet.k);
    }
    
    Hash_Map<IndexT, uint32_t> pose_to_uf;
    uint32_t uf_idx = 0;
    for (IndexT p : all_poses) pose_to_uf[p] = uf_idx++;

    UnionFind uf;
    uf.InitSets(static_cast<unsigned int>(all_poses.size()));

    // A. Add strict 5 degree triplets
    for (size_t i = 0; i < all_triplets.size(); ++i) {
      if (all_triplets[i].err <= max_angular_error) {
        uf.Union(pose_to_uf[all_triplets[i].triplet.i], pose_to_uf[all_triplets[i].triplet.j]);
        uf.Union(pose_to_uf[all_triplets[i].triplet.j], pose_to_uf[all_triplets[i].triplet.k]);
        validated_triplet_indices.insert(i);
      }
    }

    // Identify nodes that are currently in the 5 degree component graph
    std::set<IndexT> poses_in_5deg;
    for (size_t idx : validated_triplet_indices) {
       poses_in_5deg.insert(all_triplets[idx].triplet.i);
       poses_in_5deg.insert(all_triplets[idx].triplet.j);
       poses_in_5deg.insert(all_triplets[idx].triplet.k);
    }

    // Precompute bounds for all components
    std::map<uint32_t, std::pair<double, double>> comp_bounds;
    for (IndexT p : all_poses) {
      if (pose_timestamps->count(p)) {
        uint32_t r = uf.Find(pose_to_uf[p]);
        double t = pose_timestamps->at(p);
        if (comp_bounds.count(r)) {
          comp_bounds[r].first = std::min(comp_bounds[r].first, t);
          comp_bounds[r].second = std::max(comp_bounds[r].second, t);
        } else {
          comp_bounds[r] = {t, t};
        }
      }
    }

    // B. Relax constraints for temporal gaps
    std::vector<size_t> candidate_indices;
    for (size_t i = 0; i < all_triplets.size(); ++i) {
      if (all_triplets[i].err > max_angular_error) candidate_indices.push_back(i);
    }
    std::sort(candidate_indices.begin(), candidate_indices.end(), 
              [&](size_t a, size_t b) { return all_triplets[a].err < all_triplets[b].err; });

    std::vector<IndexT> rescued_nodes;

    for (size_t idx : candidate_indices) {
      const auto & tr = all_triplets[idx].triplet;
      uint32_t roots[3] = {uf.Find(pose_to_uf[tr.i]), uf.Find(pose_to_uf[tr.j]), uf.Find(pose_to_uf[tr.k])};
      
      bool merges_useful = false;
      if (roots[0] != roots[1] || roots[1] != roots[2] || roots[0] != roots[2]) {
        for (int a = 0; a < 3 && !merges_useful; ++a) {
          for (int b = a + 1; b < 3 && !merges_useful; ++b) {
            if (roots[a] == roots[b]) continue;
            if (comp_bounds.count(roots[a]) && comp_bounds.count(roots[b])) {
              double dist = std::max(0.0, std::max(comp_bounds[roots[a]].first - comp_bounds[roots[b]].second, 
                                                  comp_bounds[roots[b]].first - comp_bounds[roots[a]].second));
              if (dist <= 5.0) merges_useful = true;
            }
          }
        }
      }

      if (merges_useful) {
        validated_triplet_indices.insert(idx);
        uf.Union(roots[0], roots[1]);
        uint32_t new_root = uf.Find(roots[0]);
        uf.Union(new_root, roots[2]);
        new_root = uf.Find(new_root);

        // Update bounds for the new merged component
        std::pair<double, double> merged_b(1e18, -1e18);
        for (int i=0; i<3; ++i) {
          if (comp_bounds.count(roots[i])) {
            merged_b.first = std::min(merged_b.first, comp_bounds[roots[i]].first);
            merged_b.second = std::max(merged_b.second, comp_bounds[roots[i]].second);
          }
        }
        comp_bounds[new_root] = merged_b;
      }
    }

    // Log rescued nodes
    std::set<IndexT> poses_in_rescued;
    for (size_t idx : validated_triplet_indices) {
       poses_in_rescued.insert(all_triplets[idx].triplet.i);
       poses_in_rescued.insert(all_triplets[idx].triplet.j);
       poses_in_rescued.insert(all_triplets[idx].triplet.k);
    }
    
    for (IndexT p : poses_in_rescued) {
      if (poses_in_5deg.count(p) == 0) {
        rescued_nodes.push_back(p);
        if (pose_img_names && pose_img_names->count(p)) {
          OPENMVG_LOG_INFO << "RESCUED critical node: " << p << " (" << pose_img_names->at(p) << ")";
        } else {
          OPENMVG_LOG_INFO << "RESCUED critical node: " << p;
        }
      }
    }
    if (!rescued_nodes.empty()) {
      OPENMVG_LOG_INFO << "Summary: Recovered " << rescued_nodes.size() << " nodes via rotation relaxation.";
    }
  }
  else
  {
    // Fallback to standard 5 degree rejection for all triplets
    for (size_t i = 0; i < all_triplets.size(); ++i) {
      if (all_triplets[i].err <= max_angular_error) {
        validated_triplet_indices.insert(i);
      }
    }
  }

  //--
  // 3. FINAL VALIDATION
  //--

  for (size_t idx : validated_triplet_indices)
  {
    const graph::Triplet & triplet = all_triplets[idx].triplet;
    vec_triplets_validated.push_back(triplet);

    const IndexT I = triplet.i, J = triplet.j , K = triplet.k;
    const Pair ij(I,J), ji(J,I), jk(J,K), kj(K,J), ki(K,I), ik(I,K);

    if (map_relatives.count(ij)) map_relatives_validated[ij] = map_relatives.at(ij);
    else if (map_relatives.count(ji)) map_relatives_validated[ji] = map_relatives.at(ji);

    if (map_relatives.count(jk)) map_relatives_validated[jk] = map_relatives.at(jk);
    else if (map_relatives.count(kj)) map_relatives_validated[kj] = map_relatives.at(kj);

    if (map_relatives.count(ki)) map_relatives_validated[ki] = map_relatives.at(ki);
    else if (map_relatives.count(ik)) map_relatives_validated[ik] = map_relatives.at(ik);
  }

  map_relatives = std::move(map_relatives_validated);

  // update to keep only useful triplets
  relativeRotations.clear();
  relativeRotations.reserve(map_relatives.size());
  std::transform(map_relatives.cbegin(), map_relatives.cend(), std::back_inserter(relativeRotations), stl::RetrieveValue());
  std::transform(map_relatives.cbegin(), map_relatives.cend(), std::inserter(used_pairs, used_pairs.begin()), stl::RetrieveKey());

  // Display statistics about rotation triplets error:
  std::ostringstream os;
  os << "Statistics about rotation triplets:\n";
  minMaxMeanMedian<float>(vec_errToIdentityPerTriplet.cbegin(), vec_errToIdentityPerTriplet.cend(), os);

  std::sort(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  if (!vec_errToIdentityPerTriplet.empty())
  {
    Histogram<float> histo(0.0f, *max_element(vec_errToIdentityPerTriplet.cbegin(), vec_errToIdentityPerTriplet.cend()), 20);
    histo.Add(vec_errToIdentityPerTriplet.cbegin(), vec_errToIdentityPerTriplet.cend());
    os << histo.ToString() << "\n";
  }

  {
    os << "\nTriplets filtering based on unit cycle rotation composition error:"
      << "\n#Triplets before: " << triplets_before_count
      << "\n#Triplets after: " << vec_triplets_validated.size();
    OPENMVG_LOG_INFO << os.str();
  }

  vec_triplets = std::move(vec_triplets_validated);

  const size_t edges_end_count = relativeRotations.size();
  OPENMVG_LOG_INFO << "\n #Edges removed by triplet inference: " << edges_start_count - edges_end_count;
}

} // namespace sfm
} // namespace openMVG
