// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2026

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_SHARDED_PAIRWISE_MATCHES_HPP
#define OPENMVG_MATCHING_SHARDED_PAIRWISE_MATCHES_HPP

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "openMVG/matching/indMatch_utils.hpp"

namespace openMVG {
namespace matching {

class ShardedPairWiseMatchesWriter : public PairWiseMatchesContainer
{
public:
  explicit ShardedPairWiseMatchesWriter(
    const std::string & manifest_filename,
    const ShardedMatchFileOptions & options = {});
  ~ShardedPairWiseMatchesWriter() override;

  bool IsOpen() const;
  bool Finalize();
  std::size_t pair_count() const;
  void SetPairWrittenCallback(const std::function<void(const Pair &)> & callback);

  void insert(std::pair<Pair, IndMatches> && pairWiseMatches) override;

private:
  struct Impl;
  Impl * impl_;
};

}  // namespace matching
}  // namespace openMVG

#endif // OPENMVG_MATCHING_SHARDED_PAIRWISE_MATCHES_HPP
