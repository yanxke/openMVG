// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2019 Pierre MOULON, Romuald PERROT

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/exif/exif_IO_EasyExif.hpp"
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
#include "openMVG/matching_image_collection/E_ACRobust_WithPriors.hpp"
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

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
/// - Compute geometric coherent feature matches (robust model estimation from putative matches)
///   with motion prior constraints integrated into RANSAC
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

  std::string  sGeometricModel   = "e"; // Default to essential matrix for priors
  bool         bForce            = false;
  bool         bGuided_matching  = false;
  int          imax_iteration    = 2048;
  unsigned int ui_max_cache_size = 0;

  // Motion prior parameters
  double yaw_tol = 15.0;     // degrees
  double pitch_tol = 10.0;   // degrees
  double roll_tol = 10.0;    // degrees
  double alt_tol = 0.2;      // translation ty
  double prior_weight = 1.0; // weight of prior penalty
  double heading_max = 140.0; // max allowed heading difference
  int    spot_sample_size = 15; // number of points to spot-check

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
  // Motion prior options
  cmd.add( make_option( 'y', yaw_tol, "yaw_tol" ) );
  cmd.add( make_option( 'x', pitch_tol, "pitch_tol" ) );
  cmd.add( make_option( 'z', roll_tol, "roll_tol" ) );
  cmd.add( make_option( 'a', alt_tol, "alt_tol" ) );
  cmd.add( make_option( 'w', prior_weight, "prior_weight" ) );
  cmd.add( make_option( 'H', heading_max, "heading_max" ) );
  cmd.add( make_option( 'S', spot_sample_size, "spot_sample_size" ) );

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
                     << "   f: fundamental matrix,\n"
                     << "   e: (default) essential matrix with motion priors,\n"
                     << "   h: homography matrix.\n"
                     << "   a: essential matrix with an angular parametrization,\n"
                     << "   u: upright essential matrix with an angular parametrization,\n"
                     << "   o: orthographic essential matrix.\n"
                     << "[-r|--guided_matching]  Use the found model to improve the pairwise correspondences.\n"
                     << "[-c|--cache_size]\n"
                     << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
                     << "  If not used, all regions will be load in memory.\n"
                     << "\n[Motion Prior Options] (only apply with -g e)\n"
                     << "[-y|--yaw_tol]          Yaw tolerance in degrees (default: 15.0)\n"
                     << "[-x|--pitch_tol]        Pitch tolerance in degrees (default: 10.0)\n"
                     << "[-z|--roll_tol]         Roll tolerance in degrees (default: 10.0)\n"
                     << "[-a|--alt_tol]          Altitude (ty) tolerance (default: 0.2)\n"
                     << "[-w|--prior_weight]     Weight of prior penalty in RANSAC (default: 1.0)\n"
                     << "[-H|--heading_max]      Max allowed heading difference (default: 140.0, 180.0 to disable)\n"
                     << "[-S|--spot_sample_size] Number of points to spot-check (default: 15, 0 to disable)";

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
                   << "--cache_size         " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << "\n"
                   << "Motion priors:       "
                   << "\n"
                   << "--yaw_tol            " << yaw_tol << "\n"
                   << "--pitch_tol          " << pitch_tol << "\n"
                   << "--roll_tol           " << roll_tol << "\n"
                    << "--alt_tol            " << alt_tol << "\n"
                    << "--prior_weight       " << prior_weight << "\n"
                    << "--heading_max        " << heading_max << "\n"
                    << "--spot_sample_size   " << spot_sample_size;

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

  EGeometricModel eGeometricModelToCompute = ESSENTIAL_MATRIX; // Default for priors
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
  // Cache EXIF headings for motion priors
  //---------------------------------------
  std::map<IndexT, double> map_headings;
  if (eGeometricModelToCompute == ESSENTIAL_MATRIX)
  {
    OPENMVG_LOG_INFO << "Caching GPS headings for motion priors...";
    for (const auto & view_ptr : sfm_data.GetViews())
    {
      const std::string image_path = stlplus::folder_append_separator(sfm_data.s_root_path) + view_ptr.second->s_Img_path;
      std::unique_ptr<openMVG::exif::Exif_IO> exifIO(new openMVG::exif::Exif_IO_EasyExif(image_path));
      double heading;
      if (exifIO->GPSImageDirection(&heading))
      {
        map_headings[view_ptr.first] = heading;
      }
    }
    OPENMVG_LOG_INFO << "Loaded GPS headings for " << map_headings.size() << " / " << sfm_data.GetViews().size() << " images.";

    // Print histogram of headings (36 buckets of 10 degrees each)
    if (!map_headings.empty())
    {
      std::vector<int> heading_histogram(36, 0);
      for (const auto & h : map_headings)
      {
        double heading_val = h.second;
        // Normalize to [0, 360)
        while (heading_val < 0) heading_val += 360.0;
        while (heading_val >= 360.0) heading_val -= 360.0;
        int bucket = static_cast<int>(heading_val / 10.0);
        if (bucket >= 0 && bucket < 36)
          heading_histogram[bucket]++;
      }
      OPENMVG_LOG_INFO << "GPS Heading histogram (10-degree buckets):";
      for (int i = 0; i < 36; ++i)
      {
        if (heading_histogram[i] > 0)
        {
          OPENMVG_LOG_INFO << "  [" << (i * 10) << "-" << ((i + 1) * 10) << "): " << heading_histogram[i];
        }
      }
    }
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
  system::LoggerProgress progress(1, {}, 1);

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

  if ( !sInputPairsFilename.empty() )
  {
    // Load input pairs
    OPENMVG_LOG_INFO << "Loading input pairs ...";
    Pair_Set input_pairs;
    loadPairs( sfm_data.GetViews().size(), sInputPairsFilename, input_pairs );

    // Filter matches with the given pairs
    OPENMVG_LOG_INFO << "Filtering matches with the given pairs.";
    map_PutativeMatches = getPairs( map_PutativeMatches, input_pairs );
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
        // Use the custom filter with motion priors
        MotionPriorConfig prior_config;
        prior_config.pitch_tol = pitch_tol;
        prior_config.roll_tol = roll_tol;
        prior_config.yaw_tol = yaw_tol;
        prior_config.alt_tol = alt_tol;
        prior_config.prior_weight = prior_weight;
        prior_config.heading_max = heading_max;
        prior_config.spot_sample_size = spot_sample_size;

        filter_ptr->Robust_model_estimation(
            GeometricFilter_EMatrix_AC_WithPriors( 4.0, imax_iteration, prior_config, &map_headings, map_PutativeMatches.size() ),
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

    // -- export Geometric View Graph statistics
    graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_GeometricMatches));

    // Print filtering summary
    OPENMVG_LOG_INFO << "\n=== Geometric Filtering Summary ===";
    OPENMVG_LOG_INFO << "Putative pairs:  " << map_PutativeMatches.size();
    OPENMVG_LOG_INFO << "Geometric pairs: " << map_GeometricMatches.size();

    if (eGeometricModelToCompute == ESSENTIAL_MATRIX)
    {
      OPENMVG_LOG_INFO << "\nMotion Prior Configuration:";
      OPENMVG_LOG_INFO << "  GPS headings loaded: " << map_headings.size() << " / " << sfm_data.GetViews().size();
      OPENMVG_LOG_INFO << "  Pitch tolerance:     " << pitch_tol << " degrees";
      OPENMVG_LOG_INFO << "  Roll tolerance:      " << roll_tol << " degrees";
      OPENMVG_LOG_INFO << "  Yaw tolerance:       " << yaw_tol << " degrees";
      OPENMVG_LOG_INFO << "  Altitude tolerance:  " << alt_tol;
      OPENMVG_LOG_INFO << "  Prior weight:        " << prior_weight;

      if (map_GeometricMatches.empty() && !map_PutativeMatches.empty())
      {
        OPENMVG_LOG_WARNING << "\n*** WARNING: All pairs were filtered out! ***";
        OPENMVG_LOG_WARNING << "This may indicate that motion priors are too strict.";
        OPENMVG_LOG_WARNING << "Suggestions:";
        OPENMVG_LOG_WARNING << "  1. Try with -w 0.0 to disable priors and verify standard filtering works";
        OPENMVG_LOG_WARNING << "  2. Increase tolerances (-x, -y, -z, -a)";
        OPENMVG_LOG_WARNING << "  3. Reduce prior weight (-w) to make priors softer (e.g., -w 0.5)";
        OPENMVG_LOG_WARNING << "  4. Check if GPS headings are valid (see histogram above)";
      }
    }

    OPENMVG_LOG_INFO << "\nTask done in (s): " << timer.elapsed();

    //-- export Adjacency matrix
    OPENMVG_LOG_INFO <<  "\n Export Adjacency Matrix of the pairwise's geometric matches";

    PairWiseMatchingToAdjacencyMatrixSVG( sfm_data.GetViews().size(),
                                          map_GeometricMatches,
                                          stlplus::create_filespec( sMatchesDirectory, "GeometricAdjacencyMatrix", "svg" ) );

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
