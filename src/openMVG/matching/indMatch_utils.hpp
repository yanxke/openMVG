// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_IND_MATCH_UTILS_HPP
#define OPENMVG_MATCHING_IND_MATCH_UTILS_HPP

#include <cstdint>
#include <functional>
#include <vector>
#include <string>

#include "openMVG/matching/indMatch.hpp"

namespace openMVG {
namespace matching {

struct ShardedMatchFileOptions
{
  std::uint64_t max_shard_size_bytes = 10ull * 1024ull * 1024ull * 1024ull;
};

using MatchReadCallback = std::function<bool(const Pair &, IndMatches &&)>;

bool Load
(
  PairWiseMatches & matches,
  const std::string & filename
);

bool Save
(
  const PairWiseMatches & matches,
  const std::string & filename
);

bool SaveJson
(
  const PairWiseMatches & matches,
  const std::string & filename
);

bool IsShardedMatchFilename(const std::string & filename);

bool SaveSharded
(
  const PairWiseMatches & matches,
  const std::string & filename,
  const ShardedMatchFileOptions & options = {}
);

bool ReadSharded
(
  const std::string & filename,
  const MatchReadCallback & callback
);

bool ReadSharded
(
  const std::string & filename,
  const Pair_Set & selected_pairs,
  const MatchReadCallback & callback
);

bool ReadShardedCamera
(
  const std::string & filename,
  const IndexT camera_id,
  const MatchReadCallback & callback
);

}  // namespace matching
}  // namespace openMVG

#endif // #define OPENMVG_MATCHING_IND_MATCH_UTILS_HPP
