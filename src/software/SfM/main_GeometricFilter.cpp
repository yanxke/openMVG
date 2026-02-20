// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2019 Pierre MOULON, Romuald PERROT

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/akaze/image_describer_akaze.hpp"
#include "openMVG/features/descriptor.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_stats.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/matching_image_collection/E_ACRobust_Angular.hpp"
#include "openMVG/matching_image_collection/Eo_Robust.hpp"
#include "openMVG/matching_image_collection/F_ACRobust.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/H_ACRobust.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/geometry/pose3.hpp"



#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <iostream>
#include <locale>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;

enum EGeometricModel
{
  FUNDAMENTAL_MATRIX       = 0,
  ESSENTIAL_MATRIX         = 1,
  HOMOGRAPHY_MATRIX        = 2,
  ESSENTIAL_MATRIX_ANGULAR = 3,
  ESSENTIAL_MATRIX_ORTHO   = 4,
  ESSENTIAL_MATRIX_UPRIGHT = 5
};

/// Write a JSON file containing the relative rotation (as quaternion) for every
/// geometrically-verified essential-matrix pair. The file is placed alongside
/// the filtered matches file and is named "matches.e.rotations.json".
///
/// The rotation is recovered by:
///   1. Lifting inlier pixel coords to bearing vectors via the camera model.
///   2. Fitting an essential matrix with the 8-point algorithm.
///   3. Decomposing E into 4 candidate poses (MotionFromEssential).
///   4. Selecting the pose with the most points passing a cheirality check.
///
/// @param map_GeometricMatches  Inlier matches keyed by (view_i, view_j).
/// @param sfm_data              SfM data providing views and intrinsics.
/// @param regions_provider      Feature regions used to look up pixel coords.
/// @param sFilteredMatchesFilename  Path to matches.e.bin (used to derive output dir).
void SavePairRotationsJson(
  const PairWiseMatches & map_GeometricMatches,
  const sfm::SfM_Data & sfm_data,
  const std::shared_ptr<sfm::Regions_Provider> & regions_provider,
  const std::string & sFilteredMatchesFilename)
{
  const std::string sRotationsJsonPath = stlplus::create_filespec(
      stlplus::folder_part(sFilteredMatchesFilename),
      "matches.e.rotations",
      "json");
  std::ofstream rot_stream(sRotationsJsonPath);
  if (!rot_stream)
  {
    OPENMVG_LOG_WARNING << "Cannot write rotations JSON to: " << sRotationsJsonPath;
    return;
  }

  rot_stream << "{\n  \"version\": 1,\n  \"pairs\": [\n";
  bool first_pair = true;

  for (const auto & match_entry : map_GeometricMatches)
  {
    const Pair & pairIdx = match_entry.first;
    const IndMatches & inlier_matches = match_entry.second;
    if (inlier_matches.size() < 5) continue;

    // Retrieve views and intrinsics
    auto vi = sfm_data.GetViews().find(pairIdx.first);
    auto vj = sfm_data.GetViews().find(pairIdx.second);
    if (vi == sfm_data.GetViews().end() || vj == sfm_data.GetViews().end()) continue;

    auto ii = sfm_data.GetIntrinsics().find(vi->second->id_intrinsic);
    auto ij = sfm_data.GetIntrinsics().find(vj->second->id_intrinsic);
    if (ii == sfm_data.GetIntrinsics().end() || ij == sfm_data.GetIntrinsics().end()) continue;

    const cameras::IntrinsicBase * cam_I = ii->second.get();
    const cameras::IntrinsicBase * cam_J = ij->second.get();
    if (!cam_I || !cam_J) continue;

    // Get pixel coordinates of inlier matches
    Mat2X xI, xJ;
    MatchesPairToMat(pairIdx, inlier_matches, &sfm_data, regions_provider, xI, xJ);
    if (xI.cols() < 5) continue;

    // Lift to bearing vectors
    const Mat3X bI = (*cam_I)(xI);
    const Mat3X bJ = (*cam_J)(xJ);

    // Fit essential matrix via 8-point algorithm (inliers are already geometrically valid)
    std::vector<Mat3> E_vec;
    openMVG::EightPointRelativePoseSolver::Solve(bI, bJ, &E_vec);
    if (E_vec.empty()) continue;

    // Decompose into up to 4 candidate poses
    std::vector<geometry::Pose3> poses;
    openMVG::MotionFromEssential(E_vec[0], &poses);
    if (poses.empty()) continue;

    // Resolve cheirality ambiguity: pick pose with most points in front of both cameras
    const geometry::Pose3 pose_identity(Mat3::Identity(), Vec3::Zero());
    int best_idx = 0, best_count = 0;
    const size_t n_test = std::min<size_t>(20, static_cast<size_t>(bI.cols()));
    for (size_t k = 0; k < poses.size(); ++k)
    {
      int count = 0;
      for (size_t m = 0; m < n_test; ++m)
      {
        Vec3 X;
        if (openMVG::Triangulate2View(
              pose_identity.rotation(), pose_identity.translation(), bI.col(m),
              poses[k].rotation(),     poses[k].translation(),      bJ.col(m),
              X, openMVG::ETriangulationMethod::DEFAULT))
        {
          ++count;
        }
      }
      if (count > best_count) { best_count = count; best_idx = static_cast<int>(k); }
    }

    // Convert rotation matrix to quaternion (w, x, y, z)
    const Mat3 R = poses[static_cast<size_t>(best_idx)].rotation();
    double qw, qx, qy, qz;
    const double tr = R(0,0) + R(1,1) + R(2,2);
    if (tr > 0.0)
    {
      const double s = std::sqrt(tr + 1.0) * 2.0; // s = 4*qw
      qw = 0.25 * s;
      qx = (R(2,1) - R(1,2)) / s;
      qy = (R(0,2) - R(2,0)) / s;
      qz = (R(1,0) - R(0,1)) / s;
    }
    else if (R(0,0) > R(1,1) && R(0,0) > R(2,2))
    {
      const double s = std::sqrt(1.0 + R(0,0) - R(1,1) - R(2,2)) * 2.0; // s = 4*qx
      qw = (R(2,1) - R(1,2)) / s;
      qx = 0.25 * s;
      qy = (R(0,1) + R(1,0)) / s;
      qz = (R(0,2) + R(2,0)) / s;
    }
    else if (R(1,1) > R(2,2))
    {
      const double s = std::sqrt(1.0 + R(1,1) - R(0,0) - R(2,2)) * 2.0; // s = 4*qy
      qw = (R(0,2) - R(2,0)) / s;
      qx = (R(0,1) + R(1,0)) / s;
      qy = 0.25 * s;
      qz = (R(1,2) + R(2,1)) / s;
    }
    else
    {
      const double s = std::sqrt(1.0 + R(2,2) - R(0,0) - R(1,1)) * 2.0; // s = 4*qz
      qw = (R(1,0) - R(0,1)) / s;
      qx = (R(0,2) + R(2,0)) / s;
      qy = (R(1,2) + R(2,1)) / s;
      qz = 0.25 * s;
    }

    if (!first_pair) rot_stream << ",\n";
    first_pair = false;
    rot_stream << "    {\"i\": " << pairIdx.first << ", \"j\": " << pairIdx.second
               << ", \"qw\": " << qw << ", \"qx\": " << qx << ", \"qy\": " << qy << ", \"qz\": " << qz << "}";
  }

  rot_stream << "\n  ]\n}\n";
  OPENMVG_LOG_INFO << "Saved pair rotations to: " << sRotationsJsonPath;
}

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
/// - Export computed data
int main( int argc, char** argv )
{
  CmdLine cmd;

  // The scene
  std::string sSfM_Data_Filename;
  // The input matches
  std::string sPutativeMatchesFilename;
  // The output matches
  std::string sFilteredMatchesFilename;
  // The input pairs
  std::string sInputPairsFilename;
  // The output pairs
  std::string sOutputPairsFilename;

  std::string  sGeometricModel   = "f";
  bool         bForce            = false;
  bool         bGuided_matching  = false;
  int          imax_iteration    = 2048;
  unsigned int ui_max_cache_size = 0;

  //required
  cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  cmd.add( make_option( 'o', sFilteredMatchesFilename, "output_file" ) );
  cmd.add( make_option( 'm', sPutativeMatchesFilename, "matches" ) );
  // Options
  cmd.add( make_option( 'p', sInputPairsFilename, "input_pairs" ) );
  cmd.add( make_option( 's', sOutputPairsFilename, "output_pairs" ) );
  cmd.add( make_option( 'g', sGeometricModel, "geometric_model" ) );
  cmd.add( make_option( 'f', bForce, "force" ) );
  cmd.add( make_option( 'r', bGuided_matching, "guided_matching" ) );
  cmd.add( make_option( 'I', imax_iteration, "max_iteration" ) );
  cmd.add( make_option( 'c', ui_max_cache_size, "cache_size" ) );

  try
  {
    if ( argc == 1 )
      throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
                     << "[-i|--input_file]       A SfM_Data file\n"
                     << "[-m|--matches]          (Input) matches filename\n"
                     << "[-o|--output_file]      (Output) filtered matches filename\n"
                     << "\n[Optional]\n"
                     << "[-p|--input_pairs]      (Input) pairs filename\n"
                     << "[-s|--output_pairs]     (Output) filtered pairs filename\n"
                     << "[-f|--force]            Force to recompute data\n"
                     << "[-g|--geometric_model]\n"
                     << "  (pairwise correspondences filtering thanks to robust model estimation):\n"
                     << "   f: (default) fundamental matrix,\n"
                     << "   e: essential matrix,\n"
                     << "   h: homography matrix.\n"
                     << "   a: essential matrix with an angular parametrization,\n"
                     << "   u: upright essential matrix with an angular parametrization,\n"
                     << "   o: orthographic essential matrix.\n"
                     << "[-r|--guided_matching]  Use the found model to improve the pairwise correspondences.\n"
                     << "[-c|--cache_size]\n"
                     << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
                     << "  If not used, all regions will be load in memory.";

    OPENMVG_LOG_INFO << s;
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO << " You called : "
                   << "\n"
                   << argv[0] << "\n"
                   << "--input_file:        " << sSfM_Data_Filename << "\n"
                   << "--matches:           " << sPutativeMatchesFilename << "\n"
                   << "--output_file:       " << sFilteredMatchesFilename << "\n"
                   << "Optional parameters: "
                   << "\n"
                   << "--input_pairs        " << sInputPairsFilename << "\n"
                   << "--output_pairs       " << sOutputPairsFilename << "\n"
                   << "--force              " << (bForce ? "true" : "false") << "\n"
                   << "--geometric_model    " << sGeometricModel << "\n"
                   << "--guided_matching    " << bGuided_matching << "\n"
                   << "--cache_size         " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size));

  if ( sFilteredMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid output file";
    return EXIT_FAILURE;
  }
  if ( sSfM_Data_Filename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid SfM file";
    return EXIT_FAILURE;
  }
  if ( sPutativeMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "It is an invalid putative matche file";
    return EXIT_FAILURE;
  }

  if (!bForce && stlplus::file_exists(sFilteredMatchesFilename))
  {
    OPENMVG_LOG_INFO << "Filtered matches already exist, skipping: "
                     << sFilteredMatchesFilename;
    return EXIT_SUCCESS;
  }

  const std::string sMatchesDirectory = stlplus::folder_part( sPutativeMatchesFilename );

  EGeometricModel eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
  switch ( std::tolower(sGeometricModel[ 0 ], std::locale()) )
  {
    case 'f':
      eGeometricModelToCompute = FUNDAMENTAL_MATRIX;
      break;
    case 'e':
      eGeometricModelToCompute = ESSENTIAL_MATRIX;
      break;
    case 'h':
      eGeometricModelToCompute = HOMOGRAPHY_MATRIX;
      break;
    case 'a':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ANGULAR;
      break;
    case 'u':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_UPRIGHT;
      break;
    case 'o':
      eGeometricModelToCompute = ESSENTIAL_MATRIX_ORTHO;
      break;
    default:
      OPENMVG_LOG_ERROR << "Unknown geometric model";
      return EXIT_FAILURE;
  }

  // -----------------------------
  // - Load SfM_Data Views & intrinsics data
  // a. Load putative descriptor matches
  // [a.1] Filter matches with input pairs
  // b. Geometric filtering of putative matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfM_Data_Filename, ESfM_Data( VIEWS | INTRINSICS ) ) )
  {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \"" << sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  // Consider that the image_describer.json is inside the matches directory (which is bellow the sfm_data.bin)
  const std::string        sImage_describer = stlplus::create_filespec( sMatchesDirectory, "image_describer.json" );
  std::unique_ptr<Regions> regions_type     = Init_region_type_from_file( sImage_describer );
  if ( !regions_type )
  {
    OPENMVG_LOG_ERROR << "Invalid: " << sImage_describer << " regions type file.";
    return EXIT_FAILURE;
  }

  //---------------------------------------
  // a. Compute putative descriptor matches
  //    - Descriptor matching (according user method choice)
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------

  // Load the corresponding view regions
  std::shared_ptr<Regions_Provider> regions_provider;
  if ( ui_max_cache_size == 0 )
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>( ui_max_cache_size );
  }

  // Show the progress on the command line:
  system::LoggerProgress progress(1, {}, 5);

  if ( !regions_provider->load( sfm_data, sMatchesDirectory, regions_type, &progress ) )
  {
    OPENMVG_LOG_ERROR << "Invalid regions.";
    return EXIT_FAILURE;
  }

  PairWiseMatches map_PutativeMatches;
  //---------------------------------------
  // A. Load initial matches
  //---------------------------------------
  if ( !Load( map_PutativeMatches, sPutativeMatchesFilename ) )
  {
    OPENMVG_LOG_ERROR << "Failed to load the initial matches file.";
    return EXIT_FAILURE;
  }
  OPENMVG_LOG_INFO << "Loaded " << map_PutativeMatches.size() << " putative matching pairs.";

  if ( !sInputPairsFilename.empty() )
  {
    // Load input pairs
    OPENMVG_LOG_INFO << "Loading input pairs ...";
    Pair_Set input_pairs;
    loadPairs( sfm_data.GetViews().size(), sInputPairsFilename, input_pairs );

    // Filter matches with the given pairs
    map_PutativeMatches = getPairs( map_PutativeMatches, input_pairs );
    OPENMVG_LOG_INFO << "Number of putative pairs after input-pair filtering: " << map_PutativeMatches.size();
  }

  //---------------------------------------
  // b. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------

  std::unique_ptr<ImageCollectionGeometricFilter> filter_ptr(
      new ImageCollectionGeometricFilter( &sfm_data, regions_provider ) );

  if ( filter_ptr )
  {
    system::Timer timer;
    const double  d_distance_ratio = 0.6;

    PairWiseMatches map_GeometricMatches;
    switch ( eGeometricModelToCompute )
    {
      case HOMOGRAPHY_MATRIX:
      {
        const bool bGeometric_only_guided_matching = true;
        filter_ptr->Robust_model_estimation(
            GeometricFilter_HMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            bGeometric_only_guided_matching ? -1.0 : d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case FUNDAMENTAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_FMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_EMatrix_AC( 4.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();

        //-- Perform an additional check to remove pairs with poor overlap
        std::vector<PairWiseMatches::key_type> vec_toRemove;
        for ( const auto& pairwisematches_it : map_GeometricMatches )
        {
          const size_t putativePhotometricCount = map_PutativeMatches.find( pairwisematches_it.first )->second.size();
          const size_t putativeGeometricCount   = pairwisematches_it.second.size();
          const float  ratio                    = putativeGeometricCount / static_cast<float>( putativePhotometricCount );
          if ( putativeGeometricCount < 50 || ratio < .3f )
          {
            // the pair will be removed
            vec_toRemove.push_back( pairwisematches_it.first );
          }
        }
        //-- remove discarded pairs
        for ( const auto& pair_to_remove_it : vec_toRemove )
        {
          map_GeometricMatches.erase( pair_to_remove_it );
        }
      }
      break;
      case ESSENTIAL_MATRIX_ANGULAR:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular<false>(4.0, imax_iteration),
          map_PutativeMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_UPRIGHT:
      {
        filter_ptr->Robust_model_estimation(
          GeometricFilter_ESphericalMatrix_AC_Angular<true>(4.0, imax_iteration),
          map_PutativeMatches, bGuided_matching, d_distance_ratio, &progress);
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
      case ESSENTIAL_MATRIX_ORTHO:
      {
        filter_ptr->Robust_model_estimation(
            GeometricFilter_EOMatrix_RA( 2.0, imax_iteration ),
            map_PutativeMatches,
            bGuided_matching,
            d_distance_ratio,
            &progress );
        map_GeometricMatches = filter_ptr->Get_geometric_matches();
      }
      break;
    }

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
    if ( !Save( map_GeometricMatches, sFilteredMatchesFilename ) )
    {
      OPENMVG_LOG_ERROR << "Cannot save filtered matches in: " << sFilteredMatchesFilename;
      return EXIT_FAILURE;
    }

    // Also save as JSON for easier parsing
    const std::string sJsonFilename = stlplus::create_filespec(
      stlplus::folder_part(sFilteredMatchesFilename),
      stlplus::basename_part(sFilteredMatchesFilename),
      "json");
    if ( !SaveJson( map_GeometricMatches, sJsonFilename ) )
    {
      OPENMVG_LOG_WARNING << "Cannot save JSON matches to: " << sJsonFilename;
    }
    else
    {
      OPENMVG_LOG_INFO << "Saved JSON matches to: " << sJsonFilename;
    }

    // Write per-pair relative rotations for ESSENTIAL_MATRIX model types.
    if (eGeometricModelToCompute == ESSENTIAL_MATRIX ||
        eGeometricModelToCompute == ESSENTIAL_MATRIX_UPRIGHT ||
        eGeometricModelToCompute == ESSENTIAL_MATRIX_ANGULAR)
    {
      SavePairRotationsJson(map_GeometricMatches, sfm_data, regions_provider, sFilteredMatchesFilename);
    }


    // -- export Geometric View Graph statistics
    graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_GeometricMatches));

    OPENMVG_LOG_INFO << "Task done in (s): " << timer.elapsed();

    const Pair_Set outputPairs = getPairs( map_GeometricMatches );

    //-- export view pair graph once geometric filter have been done
    {
      std::set<IndexT> set_ViewIds;
      std::transform( sfm_data.GetViews().begin(), sfm_data.GetViews().end(), std::inserter( set_ViewIds, set_ViewIds.begin() ), stl::RetrieveKey() );
      graph::indexedGraph putativeGraph( set_ViewIds, outputPairs );
      graph::exportToGraphvizData(
          stlplus::create_filespec( sMatchesDirectory, "geometric_matches" ),
          putativeGraph );
    }

    // Write pairs
    if ( !sOutputPairsFilename.empty() )
    {
      OPENMVG_LOG_INFO << "Saving pairs to: " << sOutputPairsFilename;
      if ( !savePairs( sOutputPairsFilename, outputPairs ) )
      {
        OPENMVG_LOG_ERROR << "Failed to write pairs file";
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}
