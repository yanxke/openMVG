// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2019 Romuald PERROT

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching_image_collection/Pair_Builder.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @brief Current list of available pair mode
 *
 */
enum EPairMode
{
  PAIR_EXHAUSTIVE = 0, // Build every combination of image pairs
  PAIR_CONTIGUOUS = 1, // Only consecutive image pairs (useful for video mode)
  PAIR_COMPASS = 2     // Only pairs with compass heading within threshold (requires GPS heading in EXIF)
};

using namespace openMVG;
using namespace openMVG::sfm;

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << '\n'
            << "[-i|--input_file]         A SfM_Data file\n"
            << "[-o|--output_file]        Output file where pairs are stored\n"
            << "\n[Optional]\n"
            << "[-m|--pair_mode] mode     Pair generation mode\n"
            << "       EXHAUSTIVE:        Build all possible pairs. [default]\n"
            << "       CONTIGUOUS:        Build pairs for contiguous images (use it with --contiguous_count parameter)\n"
            << "       COMPASS:           Build pairs where compass heading difference <= threshold\n"
            << "[-c|--contiguous_count] X Number of contiguous links\n"
            << "       X: will match 0 with (1->X), ...]\n"
            << "       2: will match 0 with (1,2), 1 with (2,3), ...\n"
            << "       3: will match 0 with (1,2,3), 1 with (2,3,4), ...\n"
            << "[-h|--heading_threshold] T Maximum heading difference in degrees (default: 120.0, used with COMPASS mode)\n"
            << std::endl;
}

// This executable computes pairs of images to be matched
int main( int argc, char** argv )
{
  CmdLine cmd;

  std::string sSfMDataFilename;
  std::string sOutputPairsFilename;
  std::string sPairMode        = "EXHAUSTIVE";
  int         iContiguousCount = -1;
  double      dHeadingThreshold = 120.0;

  // Mandatory elements:
  cmd.add( make_option( 'i', sSfMDataFilename, "input_file" ) );
  cmd.add( make_option( 'o', sOutputPairsFilename, "output_file" ) );
  // Optional elements:
  cmd.add( make_option( 'm', sPairMode, "pair_mode" ) );
  cmd.add( make_option( 'c', iContiguousCount, "contiguous_count" ) );
  cmd.add( make_option( 'h', dHeadingThreshold, "heading_threshold" ) );

  try
  {
    if ( argc == 1 )
      throw std::string( "Invalid command line parameter." );
    cmd.process( argc, argv );
  }
  catch ( const std::string& s )
  {
    usage( argv[ 0 ] );
    std::cerr << "[Error] " << s << std::endl;

    return EXIT_FAILURE;
  }

  // 0. Parse parameters
  std::cout << " You called:\n"
            << argv[ 0 ] << "\n"
            << "--input_file         : " << sSfMDataFilename << "\n"
            << "--output_file        : " << sOutputPairsFilename << "\n"
            << "Optional parameters\n"
            << "--pair_mode          : " << sPairMode << "\n"
            << "--contiguous_count   : " << iContiguousCount << "\n"
            << "--heading_threshold  : " << dHeadingThreshold << "\n"
            << std::endl;

  if ( sSfMDataFilename.empty() )
  {
    usage( argv[ 0 ] );
    std::cerr << "[Error] Input file not set." << std::endl;
    exit( EXIT_FAILURE );
  }
  if ( sOutputPairsFilename.empty() )
  {
    usage( argv[ 0 ] );
    std::cerr << "[Error] Output file not set." << std::endl;
    exit( EXIT_FAILURE );
  }

  EPairMode pairMode;
  if ( sPairMode == "EXHAUSTIVE" )
  {
    pairMode = PAIR_EXHAUSTIVE;
  }
  else if ( sPairMode == "CONTIGUOUS" )
  {
    if ( iContiguousCount == -1 )
    {
      usage( argv[ 0 ] );
      std::cerr << "[Error] Contiguous pair mode selected but contiguous_count not set." << std::endl;
      exit( EXIT_FAILURE );
    }

    pairMode = PAIR_CONTIGUOUS;
  }
  else if ( sPairMode == "COMPASS" )
  {
    pairMode = PAIR_COMPASS;
  }
  else
  {
    usage( argv[ 0 ] );
    std::cerr << "[Error] Unknown pair mode: " << sPairMode << std::endl;
    exit( EXIT_FAILURE );
  }

  // 1. Load SfM data scene
  std::cout << "Loading scene.";
  SfM_Data sfm_data;
  if ( !Load( sfm_data, sSfMDataFilename, ESfM_Data( VIEWS | INTRINSICS ) ) )
  {
    std::cerr << std::endl
              << "The input SfM_Data file \"" << sSfMDataFilename << "\" cannot be read." << std::endl;
    exit( EXIT_FAILURE );
  }
  const size_t NImage = sfm_data.GetViews().size();

  // 2. Compute pairs
  std::cout << "Computing pairs." << std::endl;
  Pair_Set pairs;
  switch ( pairMode )
  {
    case PAIR_EXHAUSTIVE:
    {
      pairs = exhaustivePairs( NImage );
      break;
    }
    case PAIR_CONTIGUOUS:
    {
      pairs = contiguousWithOverlap( NImage, iContiguousCount );
      break;
    }
    case PAIR_COMPASS:
    {
      // Build pairs based on GPS compass heading
      // First, extract headings from views
      std::map<IndexT, double> headings;
      size_t views_with_heading = 0;

      for (const auto & view_pair : sfm_data.GetViews())
      {
        const IndexT view_id = view_pair.first;
        const auto view_ptr = view_pair.second;

        // Try to cast to ViewPriors to check for heading
        const ViewPriors * view_priors = dynamic_cast<const ViewPriors*>(view_ptr.get());
        if (view_priors && view_priors->b_has_heading_)
        {
          headings[view_id] = view_priors->gps_heading_;
          views_with_heading++;
        }
      }

      std::cout << "Found GPS headings for " << views_with_heading << " / " << NImage << " images." << std::endl;
      std::cout << "Heading threshold: " << dHeadingThreshold << " degrees" << std::endl;

      if (views_with_heading < 2)
      {
        std::cerr << "[Warning] Insufficient GPS headings for COMPASS mode. Falling back to EXHAUSTIVE." << std::endl;
        pairs = exhaustivePairs( NImage );
      }
      else
      {
        // Generate pairs where heading difference is within threshold
        size_t rejected = 0;
        for (size_t i = 0; i < NImage; ++i)
        {
          for (size_t j = i + 1; j < NImage; ++j)
          {
            // Check if both images have headings
            auto it_i = headings.find(i);
            auto it_j = headings.find(j);

            if (it_i != headings.end() && it_j != headings.end())
            {
              // Compute heading difference (handle wrap-around at 0/360)
              double diff = std::abs(it_i->second - it_j->second);
              if (diff > 180.0)
                diff = 360.0 - diff;

              if (diff <= dHeadingThreshold)
              {
                pairs.insert( std::make_pair( i, j ) );
              }
              else
              {
                rejected++;
              }
            }
            else
            {
              // If either image lacks heading, include the pair (conservative)
              pairs.insert( std::make_pair( i, j ) );
            }
          }
        }

        const size_t total = pairs.size() + rejected;
        const double rejection_ratio = total > 0 ? (100.0 * rejected / total) : 0.0;
        std::cout << "Generated " << pairs.size() << " pairs (rejected " << rejected << " based on heading, ratio rejected: " << std::fixed << std::setprecision(1) << rejection_ratio << "%)" << std::endl;
      }
      break;
    }
    default:
    {
      std::cerr << "Unknown pair mode" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  // 3. Save pairs
  std::cout << "Saving pairs." << std::endl;
  if ( !savePairs( sOutputPairsFilename, pairs ) )
  {
    std::cerr << "Failed to save pairs to file: \"" << sOutputPairsFilename << "\"" << std::endl;
    exit( EXIT_FAILURE );
  }

  return EXIT_SUCCESS;
}
