// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2019 Pierre MOULON, Romuald Perrot.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/graph/graph.hpp"
#include "openMVG/graph/graph_stats.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/pairwiseAdjacencyDisplay.hpp"
#include "openMVG/matching/sharded_pairwise_matches.hpp"
#include "openMVG/matching_image_collection/Cascade_Hashing_Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Matcher_Regions.hpp"
#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_preemptive_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider.hpp"
#include "openMVG/sfm/pipelines/sfm_regions_provider_cache.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/stl/stl.hpp"
#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;

namespace {

constexpr std::size_t kPairChunkSize = 500000;

class ForwardingMatchContainer : public PairWiseMatchesContainer
{
public:
  ForwardingMatchContainer(
    PairWiseMatchesContainer & sink,
    std::ofstream * pair_stream,
    const int min_match_count)
    : sink_(sink),
      pair_stream_(pair_stream),
      min_match_count_(min_match_count)
  {}

  void insert(std::pair<Pair, IndMatches> && pairWiseMatches) override
  {
    if (min_match_count_ > 0 &&
        static_cast<int>(pairWiseMatches.second.size()) < min_match_count_)
    {
      return;
    }

    if (pair_stream_)
    {
      (*pair_stream_) << pairWiseMatches.first.first << " " << pairWiseMatches.first.second << "\n";
    }
    ++pair_count_;
    sink_.insert(std::move(pairWiseMatches));
  }

  std::size_t pair_count() const
  {
    return pair_count_;
  }

private:
  PairWiseMatchesContainer & sink_;
  std::ofstream * pair_stream_;
  int min_match_count_ = 0;
  std::size_t pair_count_ = 0;
};

template <typename FlushFn>
bool ForEachPairChunkFromFile(
  const size_t N,
  const std::string & sFileName,
  const FlushFn & flush_fn)
{
  std::ifstream in(sFileName);
  if (!in)
  {
    OPENMVG_LOG_ERROR
      << "loadPairs: Impossible to read the specified file: \"" << sFileName << "\".";
    return false;
  }

  Pair_Set chunk_pairs;
  std::string sValue;
  std::vector<std::string> vec_str;
  while (std::getline(in, sValue))
  {
    vec_str.clear();
    stl::split(sValue, ' ', vec_str);
    const IndexT str_size(vec_str.size());
    if (str_size < 2)
    {
      OPENMVG_LOG_ERROR << "loadPairs: Invalid input file: \"" << sFileName << "\".";
      return false;
    }
    std::stringstream oss;
    oss.clear();
    oss.str(vec_str[0]);
    IndexT I = 0, J = 0;
    oss >> I;
    for (IndexT i = 1; i < str_size; ++i)
    {
      oss.clear();
      oss.str(vec_str[i]);
      oss >> J;
      if (I > N - 1 || J > N - 1)
      {
        OPENMVG_LOG_ERROR
          << "loadPairs: Invalid input file. Image out of range. "
          << "I: " << I << " J:" << J << " N:" << N << "\n"
          << "File: \"" << sFileName << "\".";
        return false;
      }
      if (I == J)
      {
        OPENMVG_LOG_ERROR << "loadPairs: Invalid input file. Image " << I
                          << " see itself. File: \"" << sFileName << "\".";
        return false;
      }
      chunk_pairs.insert({std::min(I, J), std::max(I, J)});
      if (chunk_pairs.size() >= kPairChunkSize)
      {
        if (!flush_fn(chunk_pairs))
        {
          return false;
        }
        chunk_pairs.clear();
      }
    }
  }

  if (!chunk_pairs.empty() && !flush_fn(chunk_pairs))
  {
    return false;
  }
  return true;
}

template <typename FlushFn>
bool ForEachExhaustivePairChunk(const size_t N, const FlushFn & flush_fn)
{
  Pair_Set chunk_pairs;
  for (IndexT I = 0; I < static_cast<IndexT>(N); ++I)
  {
    for (IndexT J = I + 1; J < static_cast<IndexT>(N); ++J)
    {
      chunk_pairs.insert({I, J});
      if (chunk_pairs.size() >= kPairChunkSize)
      {
        if (!flush_fn(chunk_pairs))
        {
          return false;
        }
        chunk_pairs.clear();
      }
    }
  }

  if (!chunk_pairs.empty() && !flush_fn(chunk_pairs))
  {
    return false;
  }
  return true;
}

}  // namespace

/// Compute corresponding features between a series of views:
/// - Load view images description (regions: features & descriptors)
/// - Compute putative local feature matches (descriptors matching)
int main( int argc, char** argv )
{
  CmdLine cmd;

  std::string  sSfM_Data_Filename;
  std::string  sOutputMatchesFilename = "";
  float        fDistRatio             = 0.8f;
  std::string  sPredefinedPairList    = "";
  std::string  sNearestMatchingMethod = "AUTO";
  bool         bForce                 = false;
  unsigned int ui_max_cache_size      = 0;
  double       shard_size_gb          = 10.0;

  // Pre-emptive matching parameters
  unsigned int ui_preemptive_feature_count = 200;
  double preemptive_matching_percentage_threshold = 0.08;

  //required
  cmd.add( make_option( 'i', sSfM_Data_Filename, "input_file" ) );
  cmd.add( make_option( 'o', sOutputMatchesFilename, "output_file" ) );
  cmd.add( make_option( 'p', sPredefinedPairList, "pair_list" ) );
  // Options
  cmd.add( make_option( 'r', fDistRatio, "ratio" ) );
  cmd.add( make_option( 'n', sNearestMatchingMethod, "nearest_matching_method" ) );
  cmd.add( make_option( 'f', bForce, "force" ) );
  cmd.add( make_option( 'c', ui_max_cache_size, "cache_size" ) );
  cmd.add( make_option( 'S', shard_size_gb, "shard_size_gb" ) );
  // Pre-emptive matching
  cmd.add( make_option( 'P', ui_preemptive_feature_count, "preemptive_feature_count") );


  try
  {
    if ( argc == 1 )
      throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    OPENMVG_LOG_INFO
      << "Usage: " << argv[ 0 ] << '\n'
      << "[-i|--input_file]   A SfM_Data file\n"
      << "[-o|--output_file]  Output file where computed matches are stored\n"
      << "[-p|--pair_list]    Pairs list file\n"
      << "\n[Optional]\n"
      << "[-f|--force] Force to recompute data]\n"
      << "[-r|--ratio] Distance ratio to discard non meaningful matches\n"
      << "   0.8: (default).\n"
      << "[-n|--nearest_matching_method]\n"
      << "  AUTO: auto choice from regions type,\n"
      << "  For Scalar based regions descriptor:\n"
      << "    BRUTEFORCEL2: L2 BruteForce matching,\n"
      << "    HNSWL2: L2 Approximate Matching with Hierarchical Navigable Small World graphs,\n"
      << "    HNSWL1: L1 Approximate Matching with Hierarchical Navigable Small World graphs\n"
      << "      tailored for quantized and histogram based descriptors (e.g uint8 RootSIFT)\n"
      << "    ANNL2: L2 Approximate Nearest Neighbor matching,\n"
      << "    CASCADEHASHINGL2: L2 Cascade Hashing matching.\n"
      << "    FASTCASCADEHASHINGL2: (default)\n"
      << "      L2 Cascade Hashing with precomputed hashed regions\n"
      << "     (faster than CASCADEHASHINGL2 but use more memory).\n"
      << "  For Binary based descriptor:\n"
      << "    BRUTEFORCEHAMMING: BruteForce Hamming matching,\n"
      << "    HNSWHAMMING: Hamming Approximate Matching with Hierarchical Navigable Small World graphs\n"
      << "[-c|--cache_size]\n"
      << "  Use a regions cache (only cache_size regions will be stored in memory)\n"
      << "  If not used, all regions will be load in memory."
      << "\n[-S|--shard_size_gb]\n"
      << "  Approximate maximum size for each sharded match file when --output_file uses .mshm\n"
      << "\n[Pre-emptive matching:]\n"
      << "[-P|--preemptive_feature_count] <NUMBER> Number of feature used for pre-emptive matching";

    OPENMVG_LOG_INFO << s;
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO << " You called : "
            << "\n"
            << argv[ 0 ] << "\n"
            << "--input_file " << sSfM_Data_Filename << "\n"
            << "--output_file " << sOutputMatchesFilename << "\n"
            << "--pair_list " << sPredefinedPairList << "\n"
            << "Optional parameters:"
            << "\n"
            << "--force " << bForce << "\n"
            << "--ratio " << fDistRatio << "\n"
            << "--nearest_matching_method " << sNearestMatchingMethod << "\n"
            << "--cache_size " << ((ui_max_cache_size == 0) ? "unlimited" : std::to_string(ui_max_cache_size)) << "\n"
            << "--shard_size_gb " << shard_size_gb << "\n"
            << "--preemptive_feature_used/count " << cmd.used('P') << " / " << ui_preemptive_feature_count;
  if (cmd.used('P'))
  {
    OPENMVG_LOG_INFO << "--preemptive_feature_count " << ui_preemptive_feature_count;
  }

  if ( sOutputMatchesFilename.empty() )
  {
    OPENMVG_LOG_ERROR << "No output file set.";
    return EXIT_FAILURE;
  }

  // If the matches already exists, skip recomputation and loading.
  if ( !bForce && ( stlplus::file_exists( sOutputMatchesFilename ) ) )
  {
    OPENMVG_LOG_INFO
      << "\t PREVIOUS RESULTS LOADED; skipping compute and load.";
    return EXIT_SUCCESS;
  }

  // -----------------------------
  // . Load SfM_Data Views & intrinsics data
  // . Compute putative descriptor matches
  // + Export some statistics
  // -----------------------------

  //---------------------------------------
  // Read SfM Scene (image view & intrinsics data)
  //---------------------------------------
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
    OPENMVG_LOG_ERROR << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read.";
    return EXIT_FAILURE;
  }
  const std::string sMatchesDirectory = stlplus::folder_part( sOutputMatchesFilename );

  //---------------------------------------
  // Load SfM Scene regions
  //---------------------------------------
  // Init the regions_type from the image describer file (used for image regions extraction)
  using namespace openMVG::features;
  const std::string sImage_describer = stlplus::create_filespec(sMatchesDirectory, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if (!regions_type)
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
  if (ui_max_cache_size == 0)
  {
    // Default regions provider (load & store all regions in memory)
    regions_provider = std::make_shared<Regions_Provider>();
  }
  else
  {
    // Cached regions provider (load & store regions on demand)
    regions_provider = std::make_shared<Regions_Provider_Cache>(ui_max_cache_size);
  }
  // If we use pre-emptive matching, we load less regions:
  if (ui_preemptive_feature_count > 0 && cmd.used('P'))
  {
    regions_provider = std::make_shared<Preemptive_Regions_Provider>(ui_preemptive_feature_count);
  }

  // Show region loading progress in coarse 10% steps to reduce log noise.
  system::LoggerProgress regions_load_progress(1, {}, 10);

  if (!regions_provider->load(sfm_data, sMatchesDirectory, regions_type, &regions_load_progress)) {
    OPENMVG_LOG_ERROR << "Cannot load view regions from: " << sMatchesDirectory << ".";
    return EXIT_FAILURE;
  }

  // Keep detailed progress for the matching stage.
  system::LoggerProgress progress(1, {}, 1);

  PairWiseMatches map_PutativeMatches;
  const bool use_sharded_output = IsShardedMatchFilename(sOutputMatchesFilename);

  // Build some alias from SfM_Data Views data:
  // - List views as a vector of filenames & image sizes
  std::vector<std::string>               vec_fileNames;
  std::vector<std::pair<size_t, size_t>> vec_imagesSize;
  {
    vec_fileNames.reserve(sfm_data.GetViews().size());
    vec_imagesSize.reserve(sfm_data.GetViews().size());
    for (const auto view_it : sfm_data.GetViews())
    {
      const View * v = view_it.second.get();
      vec_fileNames.emplace_back(stlplus::create_filespec(sfm_data.s_root_path,
          v->s_Img_path));
      vec_imagesSize.emplace_back(v->ui_width, v->ui_height);
    }
  }

  OPENMVG_LOG_INFO << " - PUTATIVE MATCHES - ";
  // Compute the putative matches
  {
    // Allocate the right Matcher according the Matching requested method
    std::unique_ptr<Matcher> collectionMatcher;
    if ( sNearestMatchingMethod == "AUTO" )
    {
      if ( regions_type->IsScalar() )
      {
        OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
        collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
      }
      else
      if (regions_type->IsBinary())
      {
        OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
        collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
      }
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEL2")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_L2));
    }
    else
    if (sNearestMatchingMethod == "BRUTEFORCEHAMMING")
    {
      OPENMVG_LOG_INFO << "Using BRUTE_FORCE_HAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, BRUTE_FORCE_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "HNSWL2")
    {
      OPENMVG_LOG_INFO << "Using HNSWL2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L2));
    }
    if (sNearestMatchingMethod == "HNSWL1")
    {
      OPENMVG_LOG_INFO << "Using HNSWL1 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_L1));
    }
    else
    if (sNearestMatchingMethod == "HNSWHAMMING")
    {
      OPENMVG_LOG_INFO << "Using HNSWHAMMING matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, HNSW_HAMMING));
    }
    else
    if (sNearestMatchingMethod == "ANNL2")
    {
      OPENMVG_LOG_INFO << "Using ANN_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, ANN_L2));
    }
    else
    if (sNearestMatchingMethod == "CASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Matcher_Regions(fDistRatio, CASCADE_HASHING_L2));
    }
    else
    if (sNearestMatchingMethod == "FASTCASCADEHASHINGL2")
    {
      OPENMVG_LOG_INFO << "Using FAST_CASCADE_HASHING_L2 matcher";
      collectionMatcher.reset(new Cascade_Hashing_Matcher_Regions(fDistRatio));
    }
    if (!collectionMatcher)
    {
      OPENMVG_LOG_ERROR << "Invalid Nearest Neighbor method: " << sNearestMatchingMethod;
      return EXIT_FAILURE;
    }
    // Perform the matching
    system::Timer timer;
    {
      // Photometric matching of putative pairs
      const std::string sOutputPairFilename =
        stlplus::create_filespec( sMatchesDirectory, "preemptive_pairs", "txt" );
      std::ofstream pair_stream(sOutputPairFilename.c_str());
      if (!pair_stream)
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches pairs in: "
          << sOutputPairFilename;
        return EXIT_FAILURE;
      }

      const int match_count_threshold = cmd.used('P')
        ? static_cast<int>(preemptive_matching_percentage_threshold * ui_preemptive_feature_count)
        : 0;

      std::unique_ptr<ShardedPairWiseMatchesWriter> sharded_writer;
      std::unique_ptr<ForwardingMatchContainer> forwarding_container;
      if (use_sharded_output)
      {
        const auto max_shard_size_bytes =
          static_cast<std::uint64_t>(std::max(1.0, shard_size_gb) * 1024.0 * 1024.0 * 1024.0);
        sharded_writer.reset(new ShardedPairWiseMatchesWriter(
          sOutputMatchesFilename,
          ShardedMatchFileOptions{max_shard_size_bytes}));
        if (!sharded_writer->IsOpen())
        {
          OPENMVG_LOG_ERROR << "Cannot initialize sharded matches output: " << sOutputMatchesFilename;
          return EXIT_FAILURE;
        }
        forwarding_container.reset(new ForwardingMatchContainer(
          *sharded_writer,
          &pair_stream,
          match_count_threshold));
      }
      else
      {
        forwarding_container.reset(new ForwardingMatchContainer(
          map_PutativeMatches,
          &pair_stream,
          match_count_threshold));
      }

      if (cmd.used('P'))
      {
        OPENMVG_LOG_INFO << "Applying preemptive match filtering during matching with threshold: "
                         << match_count_threshold;
      }

      const auto match_chunk =
        [&](const Pair_Set & pair_chunk) -> bool
        {
          if (pair_chunk.empty())
          {
            return true;
          }
          OPENMVG_LOG_INFO << "Running matching on pair chunk of size: " << pair_chunk.size();
          collectionMatcher->Match( regions_provider, pair_chunk, *forwarding_container, &progress );
          return true;
        };

      if ( sPredefinedPairList.empty() )
      {
        OPENMVG_LOG_INFO << "No input pair file set. Use exhaustive match by default.";
        const size_t NImage = sfm_data.GetViews().size();
        if (!ForEachExhaustivePairChunk(NImage, match_chunk))
        {
          return EXIT_FAILURE;
        }
      }
      else if (!ForEachPairChunkFromFile(sfm_data.GetViews().size(), sPredefinedPairList, match_chunk))
      {
        OPENMVG_LOG_ERROR << "Failed to load pairs from file: \"" << sPredefinedPairList << "\"";
        return EXIT_FAILURE;
      }

      //---------------------------------------
      //-- Export putative matches & pairs
      //---------------------------------------
      OPENMVG_LOG_INFO << "Saving putative matches...";
      if (use_sharded_output)
      {
        if (!sharded_writer->Finalize())
        {
          OPENMVG_LOG_ERROR
            << "Cannot finalize sharded matches in: "
            << sOutputMatchesFilename;
          return EXIT_FAILURE;
        }
      }
      else if ( !Save( map_PutativeMatches, std::string( sOutputMatchesFilename ) ) )
      {
        OPENMVG_LOG_ERROR
          << "Cannot save computed matches in: "
          << sOutputMatchesFilename;
        return EXIT_FAILURE;
      }
      OPENMVG_LOG_INFO << "Putative matches saved.";
      OPENMVG_LOG_INFO << "Pairs saved.";
    }
    OPENMVG_LOG_INFO << "Task (Regions Matching) done in (s): " << timer.elapsed();
  }

  OPENMVG_LOG_INFO << "#Putative pairs: "
                   << (use_sharded_output ? "stored in sharded manifest" : std::to_string(map_PutativeMatches.size()));

  if (!use_sharded_output)
  {
    // -- export Putative View Graph statistics
    graph::getGraphStatistics(sfm_data.GetViews().size(), getPairs(map_PutativeMatches));

    //-- export view pair graph once putative graph matches has been computed
    {
      std::set<IndexT> set_ViewIds;
      std::transform( sfm_data.GetViews().begin(), sfm_data.GetViews().end(), std::inserter( set_ViewIds, set_ViewIds.begin() ), stl::RetrieveKey() );
      graph::indexedGraph putativeGraph( set_ViewIds, getPairs( map_PutativeMatches ) );
      graph::exportToGraphvizData(
          stlplus::create_filespec( sMatchesDirectory, "putative_matches" ),
          putativeGraph );
    }
  }
  else
  {
    OPENMVG_LOG_WARNING
      << "Skipping in-memory graph statistics and Graphviz export for sharded putative matches.";
  }

  return EXIT_SUCCESS;
}
