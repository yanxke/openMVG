// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
// Copyright (c) 2026
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/system/logger.hpp"
#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/istreamwrapper.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;

static std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    if (!item.empty()) elems.push_back(item);
  }
  return elems;
}

struct RelRotation {

  double qw, qx, qy, qz;
};

// Custom JSON match loader (since indMatch_utils only supports counts in JSON)
bool LoadMatchesJson(const std::string & filename, PairWiseMatches & matches)
{
  std::ifstream ifs(filename.c_str());
  if (!ifs.is_open()) return false;

  rapidjson::IStreamWrapper isw(ifs);
  rapidjson::Document doc;
  doc.ParseStream(isw);
  if (doc.HasParseError() || !doc.IsObject() || !doc.HasMember("pairs") || !doc["pairs"].IsArray())
    return false;

  for (const auto & entry : doc["pairs"].GetArray())
  {
    if (!entry.IsObject() || !entry.HasMember("i") || !entry.HasMember("j")) continue;
    IndexT i = entry["i"].GetUint();
    IndexT j = entry["j"].GetUint();
    Pair p(std::min(i, j), std::max(i, j));
    
    IndMatches ind_matches;
    if (entry.HasMember("matches") && entry["matches"].IsArray())
    {
      for (const auto & m : entry["matches"].GetArray())
      {
        if (m.IsArray() && m.Size() == 2)
        {
          ind_matches.emplace_back(m[0].GetUint(), m[1].GetUint());
        }
      }
    }
    matches[p] = std::move(ind_matches);
  }
  return true;
}

// Custom JSON match saver (includes full match indices)
bool SaveMatchesJsonFull(const PairWiseMatches & matches, const std::string & filename)
{
  std::ofstream ofs(filename.c_str());
  if (!ofs.is_open()) return false;

  ofs << "{\n  \"pairs\": [\n";
  bool first_pair = true;
  for (const auto & kv : matches)
  {
    if (!first_pair) ofs << ",\n";
    first_pair = false;
    ofs << "    {\"i\": " << kv.first.first << ", \"j\": " << kv.first.second 
        << ", \"num_matches\": " << kv.second.size() << ", \"matches\": [";
    bool first_m = true;
    for (const auto & m : kv.second)
    {
      if (!first_m) ofs << ",";
      first_m = false;
      ofs << "[" << m.i_ << "," << m.j_ << "]";
    }
    ofs << "]}";
  }
  ofs << "\n  ]\n}\n";
  return true;
}


bool LoadRotationsJson(const std::string & filename, std::map<Pair, RelRotation> & rotations)
{
  std::ifstream ifs(filename.c_str());
  if (!ifs.is_open()) return false;

  rapidjson::IStreamWrapper isw(ifs);
  rapidjson::Document doc;
  doc.ParseStream(isw);
  if (doc.HasParseError() || !doc.IsObject() || !doc.HasMember("pairs") || !doc["pairs"].IsArray())
    return false;

  for (const auto & entry : doc["pairs"].GetArray())
  {
    if (!entry.IsObject() || !entry.HasMember("i") || !entry.HasMember("j") || 
        !entry.HasMember("qw") || !entry.HasMember("qx") || !entry.HasMember("qy") || !entry.HasMember("qz"))
      continue;
    IndexT i = entry["i"].GetUint();
    IndexT j = entry["j"].GetUint();
    Pair p(std::min(i, j), std::max(i, j));
    rotations[p] = {entry["qw"].GetDouble(), entry["qx"].GetDouble(), entry["qy"].GetDouble(), entry["qz"].GetDouble()};
  }
  return true;
}

bool SaveRotationsJson(const std::map<Pair, RelRotation> & rotations, const std::string & filename)
{
  std::ofstream ofs(filename.c_str());
  if (!ofs.is_open()) return false;
  ofs << "{\n  \"version\": 1,\n  \"pairs\": [\n";
  bool first = true;
  for (const auto & kv : rotations)
  {
    if (!first) ofs << ",\n";
    first = false;
    ofs << "    {\"i\": " << kv.first.first << ", \"j\": " << kv.first.second 
        << ", \"qw\": " << kv.second.qw << ", \"qx\": " << kv.second.qx 
        << ", \"qy\": " << kv.second.qy << ", \"qz\": " << kv.second.qz << "}";
  }
  ofs << "\n  ]\n}\n";
  return true;
}

int main(int argc, char ** argv)
{
  std::string sMatchFiles;
  std::string sRotationFiles;
  std::string sOutputBasename;

  CmdLine cmd;
  cmd.add(make_option('m', sMatchFiles, "matches"));
  cmd.add(make_option('r', sRotationFiles, "rotations"));
  cmd.add(make_option('o', sOutputBasename, "output_basename"));

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string & s) {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << "\n"
                     << "[-m|--matches] match_file1,match_file2,...\n"
                     << "[-r|--rotations] rot_file1,rot_file2,...\n"
                     << "[-o|--output_basename] e.g. matches.e (will save as .bin, .json, .rotations.json)\n";
    return EXIT_FAILURE;
  }

  std::vector<std::string> vec_match_files = split(sMatchFiles, ',');
  std::vector<std::string> vec_rotation_files = split(sRotationFiles, ',');

  PairWiseMatches merged_matches;
  std::map<Pair, RelRotation> merged_rotations;

  // Track match counts per pair for rotation selection
  std::map<Pair, std::vector<size_t>> pair_counts_per_input; // pair -> [count1, count2, ...]
  pair_counts_per_input.clear();

  std::vector<size_t> match_counts_per_file;
  std::vector<size_t> pair_counts_per_file;

  size_t match_file_idx = 0;
  for (const auto & f : vec_match_files)
  {
    PairWiseMatches current;
    bool ok = false;
    if (stlplus::extension_part(f) == "json") ok = LoadMatchesJson(f, current);
    else ok = Load(current, f);

    if (!ok) {
      OPENMVG_LOG_WARNING << "Failed to load matches from: " << f;
      continue;
    }

    size_t current_file_matches = 0;
    pair_counts_per_file.push_back(current.size());

    for (const auto & kv : current)
    {
      const Pair & p = kv.first;
      const IndMatches & inliers = kv.second;
      current_file_matches += inliers.size();
      
      if (merged_matches.find(p) == merged_matches.end())
        merged_matches[p] = inliers;
      else
      {
        IndMatches & merged = merged_matches[p];
        std::set<IndMatch> unique(merged.begin(), merged.end());
        for (const auto & m : inliers)
          if (unique.insert(m).second) merged.push_back(m);
      }
      
      if (pair_counts_per_input[p].size() <= match_file_idx)
        pair_counts_per_input[p].resize(match_file_idx + 1, 0);
      pair_counts_per_input[p][match_file_idx] = inliers.size();
    }
    match_counts_per_file.push_back(current_file_matches);
    match_file_idx++;
  }

  // Identify pairs in later files that were not in the first one
  size_t new_pairs_in_upright = 0;
  if (vec_match_files.size() > 1)
  {
    for (const auto & kv : pair_counts_per_input)
    {
       // If it has count in index 1 (upright) but not index 0 (standard)
       if (kv.second.size() > 1 && kv.second[1] > 0 && (kv.second[0] == 0))
          new_pairs_in_upright++;
    }
  }

  size_t rot_file_idx = 0;
  for (const auto & f : vec_rotation_files)
  {
    std::map<Pair, RelRotation> current;
    if (!LoadRotationsJson(f, current)) {
      OPENMVG_LOG_WARNING << "Failed to load rotations from: " << f;
      continue;
    }

    for (const auto & kv : current)
    {
      const Pair & p = kv.first;
      // We pick the rotation from the file that had the most inliers for this pair,
      // assuming rotation file index corresponds to match file index.
      // If we have fewer rotation files than match files, we'll need to be careful.
      if (merged_rotations.find(p) == merged_rotations.end())
        merged_rotations[p] = kv.second;
      else
      {
        // Compare current with best
        // Note: this logic assumes 1-to-1 match between -m and -r lists.
        // If not, we just pick the first one we find or the one with highest count if we can map it.
        // For simplicity, let's just keep the one with most inliers if we can find the counts.
        size_t current_count = 0;
        if (pair_counts_per_input.count(p) && pair_counts_per_input[p].size() > rot_file_idx)
          current_count = pair_counts_per_input[p][rot_file_idx];
          
        // We'll need to know which file gave the *already stored* rotation.
        // Let's just store the best count seen so far.
      }
    }
    rot_file_idx++;
  }
  
  // Re-merging rotations with "best inlier" logic
  merged_rotations.clear();
  for (size_t i = 0; i < vec_rotation_files.size(); ++i)
  {
    std::map<Pair, RelRotation> current;
    LoadRotationsJson(vec_rotation_files[i], current);
    for (const auto & kv : current)
    {
      const Pair & p = kv.first;
      if (merged_rotations.find(p) == merged_rotations.end())
        merged_rotations[p] = kv.second;
      else
      {
        // Find best count
        size_t best_count = 0;
        // Search all previous rotation files to find the best count for this pair
        for (size_t j = 0; j < i; ++j) {
           if (pair_counts_per_input.count(p) && pair_counts_per_input[p].size() > j)
             best_count = std::max(best_count, pair_counts_per_input[p][j]);
        }
        size_t current_count = 0;
        if (pair_counts_per_input.count(p) && pair_counts_per_input[p].size() > i)
          current_count = pair_counts_per_input[p][i];
        
        if (current_count > best_count) merged_rotations[p] = kv.second;
      }
    }
  }

  // Save results
  const std::string out_bin = sOutputBasename + ".bin";
  const std::string out_json = sOutputBasename + ".json";
  const std::string out_rot = sOutputBasename + ".rotations.json";

  size_t total_merged_matches = 0;
  for (const auto & kv : merged_matches) total_merged_matches += kv.second.size();

  OPENMVG_LOG_INFO << "Merge Statistics:";
  for (size_t i = 0; i < vec_match_files.size(); ++i)
  {
     OPENMVG_LOG_INFO << " - Input " << i << " (" << stlplus::basename_part(vec_match_files[i]) << "): " 
                      << pair_counts_per_file[i] << " pairs, " << match_counts_per_file[i] << " matches.";
  }
  OPENMVG_LOG_INFO << " - Merged Output: " << merged_matches.size() << " pairs, " << total_merged_matches << " matches.";
  if (vec_match_files.size() > 1)
  {
    OPENMVG_LOG_INFO << " - New pairs contributed by " << stlplus::basename_part(vec_match_files[1]) << " only: " << new_pairs_in_upright;
  }

  if (Save(merged_matches, out_bin)) OPENMVG_LOG_INFO << "Saved merged matches to: " << out_bin;
  if (SaveMatchesJsonFull(merged_matches, out_json)) OPENMVG_LOG_INFO << "Saved merged JSON matches to: " << out_json;
  if (SaveRotationsJson(merged_rotations, out_rot)) OPENMVG_LOG_INFO << "Saved merged rotations to: " << out_rot;

  return EXIT_SUCCESS;
}
