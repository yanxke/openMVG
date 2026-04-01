// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2026

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/matching/sharded_pairwise_matches.hpp"

#include "openMVG/system/logger.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace openMVG {
namespace matching {

namespace {

constexpr char kManifestMagic[] = "OPENMVG_SHARDED_MATCHES";
constexpr std::uint64_t kManifestVersion = 1;
constexpr char kShardedExtension[] = "mshm";

struct ShardManifestEntry
{
  std::string filename;
  std::uint64_t size_bytes = 0;
  std::uint64_t pair_count = 0;
};

struct ManifestData
{
  std::uint64_t max_shard_size_bytes = 0;
  std::uint64_t pair_count = 0;
  std::string camera_index_filename;
  std::vector<ShardManifestEntry> shards;
};

inline std::uint64_t ToU64(const IndexT value)
{
  return static_cast<std::uint64_t>(value);
}

inline IndexT ToIndexT(const std::uint64_t value)
{
  return static_cast<IndexT>(value);
}

inline std::uint64_t EntryByteSize(const IndMatches & matches)
{
  return 3ull * sizeof(std::uint64_t) + static_cast<std::uint64_t>(matches.size()) * 2ull * sizeof(std::uint64_t);
}

inline std::string ManifestBaseStem(const std::string & manifest_filename)
{
  const std::string folder = stlplus::folder_part(manifest_filename);
  const std::string basename = stlplus::basename_part(manifest_filename);
  if (folder.empty())
  {
    return basename;
  }
  return folder + "/" + basename;
}

inline std::string CameraIndexFilename(const std::string & manifest_filename)
{
  return manifest_filename + ".cidx";
}

inline std::string ShardFilename(const std::string & manifest_filename, const std::uint32_t shard_id)
{
  std::ostringstream os;
  os << ManifestBaseStem(manifest_filename) << ".shard.";
  os.width(5);
  os.fill('0');
  os << shard_id << ".mshard";
  return os.str();
}

bool WriteManifest(const std::string & manifest_filename, const ManifestData & manifest)
{
  std::ofstream stream(manifest_filename);
  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot write sharded matches manifest: " << manifest_filename;
    return false;
  }

  stream << kManifestMagic << " " << kManifestVersion << "\n";
  stream << "max_shard_size_bytes " << manifest.max_shard_size_bytes << "\n";
  stream << "pair_count " << manifest.pair_count << "\n";
  stream << "camera_index_file " << stlplus::filename_part(manifest.camera_index_filename) << "\n";
  stream << "shard_count " << manifest.shards.size() << "\n";
  for (const ShardManifestEntry & shard : manifest.shards)
  {
    stream << "shard "
           << stlplus::filename_part(shard.filename) << " "
           << shard.size_bytes << " "
           << shard.pair_count << "\n";
  }
  return static_cast<bool>(stream);
}

bool ReadManifest(const std::string & manifest_filename, ManifestData & manifest)
{
  manifest = ManifestData{};
  std::ifstream stream(manifest_filename);
  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot open sharded matches manifest: " << manifest_filename;
    return false;
  }

  std::string magic;
  std::uint64_t version = 0;
  if (!(stream >> magic >> version) || magic != kManifestMagic || version != kManifestVersion)
  {
    OPENMVG_LOG_ERROR << "Invalid sharded matches manifest header: " << manifest_filename;
    return false;
  }

  std::string key;
  std::size_t shard_count = 0;
  while (stream >> key)
  {
    if (key == "max_shard_size_bytes")
    {
      stream >> manifest.max_shard_size_bytes;
    }
    else if (key == "pair_count")
    {
      stream >> manifest.pair_count;
    }
    else if (key == "camera_index_file")
    {
      std::string value;
      stream >> value;
      manifest.camera_index_filename = stlplus::create_filespec(
        stlplus::folder_part(manifest_filename), value);
    }
    else if (key == "shard_count")
    {
      stream >> shard_count;
      manifest.shards.reserve(shard_count);
    }
    else if (key == "shard")
    {
      std::string shard_name;
      ShardManifestEntry entry;
      stream >> shard_name >> entry.size_bytes >> entry.pair_count;
      entry.filename = stlplus::create_filespec(
        stlplus::folder_part(manifest_filename), shard_name);
      manifest.shards.push_back(std::move(entry));
    }
  }

  if (!stream.eof())
  {
    OPENMVG_LOG_ERROR << "Failed to parse sharded matches manifest: " << manifest_filename;
    return false;
  }
  return true;
}

bool WriteCameraIndex(
  const std::string & filename,
  const std::map<IndexT, std::vector<std::uint32_t>> & camera_to_shards)
{
  std::ofstream stream(filename);
  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot write camera index: " << filename;
    return false;
  }

  stream << "OPENMVG_CAMERA_SHARD_INDEX 1\n";
  stream << "camera_count " << camera_to_shards.size() << "\n";
  for (const auto & entry : camera_to_shards)
  {
    stream << entry.first << " " << entry.second.size();
    for (const std::uint32_t shard_id : entry.second)
    {
      stream << " " << shard_id;
    }
    stream << "\n";
  }
  return static_cast<bool>(stream);
}

bool ReadCameraIndex(
  const std::string & filename,
  std::map<IndexT, std::vector<std::uint32_t>> & camera_to_shards)
{
  camera_to_shards.clear();
  std::ifstream stream(filename);
  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot open camera index: " << filename;
    return false;
  }

  std::string magic;
  std::uint64_t version = 0;
  std::size_t camera_count = 0;
  if (!(stream >> magic >> version) || magic != "OPENMVG_CAMERA_SHARD_INDEX" || version != 1)
  {
    OPENMVG_LOG_ERROR << "Invalid camera index header: " << filename;
    return false;
  }
  std::string key;
  if (!(stream >> key >> camera_count) || key != "camera_count")
  {
    OPENMVG_LOG_ERROR << "Invalid camera index body: " << filename;
    return false;
  }

  for (std::size_t i = 0; i < camera_count; ++i)
  {
    IndexT camera_id = 0;
    std::size_t shard_count = 0;
    stream >> camera_id >> shard_count;
    std::vector<std::uint32_t> shard_ids;
    shard_ids.reserve(shard_count);
    for (std::size_t j = 0; j < shard_count; ++j)
    {
      std::uint32_t shard_id = 0;
      stream >> shard_id;
      shard_ids.push_back(shard_id);
    }
    camera_to_shards.emplace(camera_id, std::move(shard_ids));
  }

  return static_cast<bool>(stream) || stream.eof();
}

bool WriteEntry(std::ofstream & shard_stream, const Pair & pair, const IndMatches & matches)
{
  const std::uint64_t first = ToU64(pair.first);
  const std::uint64_t second = ToU64(pair.second);
  const std::uint64_t match_count = static_cast<std::uint64_t>(matches.size());
  shard_stream.write(reinterpret_cast<const char *>(&first), sizeof(first));
  shard_stream.write(reinterpret_cast<const char *>(&second), sizeof(second));
  shard_stream.write(reinterpret_cast<const char *>(&match_count), sizeof(match_count));
  for (const IndMatch & match : matches)
  {
    const std::uint64_t i = ToU64(match.i_);
    const std::uint64_t j = ToU64(match.j_);
    shard_stream.write(reinterpret_cast<const char *>(&i), sizeof(i));
    shard_stream.write(reinterpret_cast<const char *>(&j), sizeof(j));
  }
  return static_cast<bool>(shard_stream);
}

bool ReadWholeShard(const std::string & shard_filename, const MatchReadCallback & callback)
{
  std::ifstream stream(shard_filename, std::ios::binary);
  if (!stream)
  {
    OPENMVG_LOG_ERROR << "Cannot open shard file: " << shard_filename;
    return false;
  }

  while (true)
  {
    std::uint64_t first = 0, second = 0, match_count = 0;
    stream.read(reinterpret_cast<char *>(&first), sizeof(first));
    if (stream.eof())
    {
      break;
    }
    stream.read(reinterpret_cast<char *>(&second), sizeof(second));
    stream.read(reinterpret_cast<char *>(&match_count), sizeof(match_count));
    if (!stream)
    {
      OPENMVG_LOG_ERROR << "Corrupted shard entry in: " << shard_filename;
      return false;
    }

    IndMatches matches;
    matches.reserve(static_cast<std::size_t>(match_count));
    for (std::uint64_t idx = 0; idx < match_count; ++idx)
    {
      std::uint64_t i = 0, j = 0;
      stream.read(reinterpret_cast<char *>(&i), sizeof(i));
      stream.read(reinterpret_cast<char *>(&j), sizeof(j));
      if (!stream)
      {
        OPENMVG_LOG_ERROR << "Corrupted shard payload in: " << shard_filename;
        return false;
      }
      matches.emplace_back(ToIndexT(i), ToIndexT(j));
    }

    if (!callback({ToIndexT(first), ToIndexT(second)}, std::move(matches)))
    {
      return false;
    }
  }
  return true;
}

}  // namespace

struct ShardedPairWiseMatchesWriter::Impl
{
  explicit Impl(const std::string & manifest_filename_, const ShardedMatchFileOptions & options_)
    : manifest_filename(manifest_filename_),
      options(options_),
      camera_index_filename(CameraIndexFilename(manifest_filename_))
  {
    const std::string output_dir = stlplus::folder_part(manifest_filename);
    if (!output_dir.empty() && !stlplus::folder_exists(output_dir))
    {
      if (!stlplus::folder_create(output_dir))
      {
        OPENMVG_LOG_ERROR << "Cannot create output directory for sharded matches: " << output_dir;
        failed = true;
        return;
      }
    }

  }

  ~Impl()
  {
    if (!finalized)
    {
      Finalize();
    }
  }

  bool OpenNextShard()
  {
    if (failed)
    {
      return false;
    }
    CloseCurrentShard();
    current_shard_id = static_cast<std::uint32_t>(shards.size());
    current_shard_filename = ShardFilename(manifest_filename, current_shard_id);
    current_shard_stream.open(current_shard_filename.c_str(), std::ios::binary);
    if (!current_shard_stream)
    {
      OPENMVG_LOG_ERROR << "Cannot create shard file: " << current_shard_filename;
      failed = true;
      return false;
    }
    current_shard_bytes = 0;
    current_shard_pairs = 0;
    return true;
  }

  void CloseCurrentShard()
  {
    if (!current_shard_stream.is_open())
    {
      return;
    }
    current_shard_stream.close();
    shards.push_back({current_shard_filename, current_shard_bytes, current_shard_pairs});
    current_shard_filename.clear();
    current_shard_bytes = 0;
    current_shard_pairs = 0;
  }

  void TouchCamera(IndexT camera_id)
  {
    std::uint32_t & last_shard = last_camera_shard[camera_id];
    if (!camera_seen.count(camera_id))
    {
      camera_to_shards[camera_id].push_back(current_shard_id);
      camera_seen.insert(camera_id);
      last_shard = current_shard_id;
      return;
    }
    if (last_shard != current_shard_id)
    {
      camera_to_shards[camera_id].push_back(current_shard_id);
      last_shard = current_shard_id;
    }
  }

  void Insert(std::pair<Pair, IndMatches> && pairWiseMatches)
  {
    if (failed)
    {
      return;
    }

    const Pair & pair = pairWiseMatches.first;
    const IndMatches & matches = pairWiseMatches.second;
    const std::uint64_t entry_bytes = EntryByteSize(matches);
    if (!current_shard_stream.is_open() ||
        (current_shard_bytes > 0 && current_shard_bytes + entry_bytes > options.max_shard_size_bytes))
    {
      if (!OpenNextShard())
      {
        return;
      }
    }

    const std::uint64_t byte_offset = current_shard_bytes;
    if (!WriteEntry(current_shard_stream, pair, matches))
    {
      OPENMVG_LOG_ERROR << "Failed to write sharded matches entry to: " << current_shard_filename;
      failed = true;
      return;
    }

    current_shard_bytes += entry_bytes;
    ++current_shard_pairs;
    ++pair_count;
    TouchCamera(pair.first);
    TouchCamera(pair.second);
    if (pair_written_callback)
    {
      pair_written_callback(pair);
    }
  }

  bool Finalize()
  {
    if (finalized)
    {
      return !failed;
    }
    finalized = true;
    if (!failed)
    {
      CloseCurrentShard();
      failed = !WriteCameraIndex(camera_index_filename, camera_to_shards);
      if (!failed)
      {
        ManifestData manifest;
        manifest.max_shard_size_bytes = options.max_shard_size_bytes;
        manifest.pair_count = pair_count;
        manifest.camera_index_filename = camera_index_filename;
        manifest.shards = shards;
        failed = !WriteManifest(manifest_filename, manifest);
      }
    }
    return !failed;
  }

  std::string manifest_filename;
  ShardedMatchFileOptions options;
  std::string camera_index_filename;
  std::ofstream current_shard_stream;
  std::string current_shard_filename;
  std::vector<ShardManifestEntry> shards;
  std::map<IndexT, std::vector<std::uint32_t>> camera_to_shards;
  std::unordered_map<IndexT, std::uint32_t> last_camera_shard;
  std::unordered_set<IndexT> camera_seen;
  std::function<void(const Pair &)> pair_written_callback;
  std::uint64_t current_shard_bytes = 0;
  std::uint64_t current_shard_pairs = 0;
  std::uint64_t pair_count = 0;
  std::uint32_t current_shard_id = 0;
  bool failed = false;
  bool finalized = false;
};

ShardedPairWiseMatchesWriter::ShardedPairWiseMatchesWriter(
  const std::string & manifest_filename,
  const ShardedMatchFileOptions & options)
  : impl_(new Impl(manifest_filename, options))
{
}

ShardedPairWiseMatchesWriter::~ShardedPairWiseMatchesWriter()
{
  delete impl_;
}

bool ShardedPairWiseMatchesWriter::IsOpen() const
{
  return impl_ && !impl_->failed;
}

bool ShardedPairWiseMatchesWriter::Finalize()
{
  return impl_ && impl_->Finalize();
}

std::size_t ShardedPairWiseMatchesWriter::pair_count() const
{
  return impl_ ? static_cast<std::size_t>(impl_->pair_count) : 0;
}

void ShardedPairWiseMatchesWriter::SetPairWrittenCallback(const std::function<void(const Pair &)> & callback)
{
  if (impl_)
  {
    impl_->pair_written_callback = callback;
  }
}

void ShardedPairWiseMatchesWriter::insert(std::pair<Pair, IndMatches> && pairWiseMatches)
{
  if (impl_)
  {
    impl_->Insert(std::move(pairWiseMatches));
  }
}

bool IsShardedMatchFilename(const std::string & filename)
{
  return stlplus::extension_part(filename) == kShardedExtension;
}

bool SaveSharded(
  const PairWiseMatches & matches,
  const std::string & filename,
  const ShardedMatchFileOptions & options)
{
  ShardedPairWiseMatchesWriter writer(filename, options);
  if (!writer.IsOpen())
  {
    return false;
  }
  for (const auto & match_entry : matches)
  {
    writer.insert({match_entry.first, match_entry.second});
  }
  return writer.Finalize();
}

bool ReadSharded(const std::string & filename, const MatchReadCallback & callback)
{
  ManifestData manifest;
  if (!ReadManifest(filename, manifest))
  {
    return false;
  }

  for (const ShardManifestEntry & shard : manifest.shards)
  {
    if (!ReadWholeShard(shard.filename, callback))
    {
      return false;
    }
  }
  return true;
}

bool ReadSharded(
  const std::string & filename,
  const Pair_Set & selected_pairs,
  const MatchReadCallback & callback)
{
  ManifestData manifest;
  if (!ReadManifest(filename, manifest))
  {
    return false;
  }
  if (selected_pairs.empty())
  {
    return ReadSharded(filename, callback);
  }

  std::map<IndexT, std::vector<std::uint32_t>> camera_to_shards;
  if (!ReadCameraIndex(manifest.camera_index_filename, camera_to_shards))
  {
    return false;
  }

  std::unordered_set<std::uint32_t> shard_ids;
  for (const Pair & pair : selected_pairs)
  {
    const auto it_first = camera_to_shards.find(pair.first);
    if (it_first != camera_to_shards.end())
    {
      shard_ids.insert(it_first->second.begin(), it_first->second.end());
    }
    const auto it_second = camera_to_shards.find(pair.second);
    if (it_second != camera_to_shards.end())
    {
      shard_ids.insert(it_second->second.begin(), it_second->second.end());
    }
  }

  std::vector<std::uint32_t> ordered_shard_ids(shard_ids.begin(), shard_ids.end());
  std::sort(ordered_shard_ids.begin(), ordered_shard_ids.end());
  for (const std::uint32_t shard_id : ordered_shard_ids)
  {
    if (shard_id >= manifest.shards.size())
    {
      OPENMVG_LOG_ERROR << "Invalid shard id in camera index: " << shard_id;
      return false;
    }

    const auto pair_filter =
      [&selected_pairs, &callback](const Pair & pair, IndMatches && matches) -> bool
      {
        if (selected_pairs.count(pair))
        {
          return callback(pair, std::move(matches));
        }
        return true;
      };

    if (!ReadWholeShard(manifest.shards[shard_id].filename, pair_filter))
    {
      return false;
    }
  }
  return true;
}

bool ReadShardedCamera(
  const std::string & filename,
  const IndexT camera_id,
  const MatchReadCallback & callback)
{
  ManifestData manifest;
  if (!ReadManifest(filename, manifest))
  {
    return false;
  }

  std::map<IndexT, std::vector<std::uint32_t>> camera_to_shards;
  if (!ReadCameraIndex(manifest.camera_index_filename, camera_to_shards))
  {
    return false;
  }

  const auto it = camera_to_shards.find(camera_id);
  if (it == camera_to_shards.end())
  {
    return true;
  }

  for (const std::uint32_t shard_id : it->second)
  {
    if (shard_id >= manifest.shards.size())
    {
      OPENMVG_LOG_ERROR << "Invalid shard id in camera index: " << shard_id;
      return false;
    }

    const auto camera_filter = [&callback, camera_id](const Pair & pair, IndMatches && matches) -> bool
    {
      if (pair.first == camera_id || pair.second == camera_id)
      {
        return callback(pair, std::move(matches));
      }
      return true;
    };

    if (!ReadWholeShard(manifest.shards[shard_id].filename, camera_filter))
    {
      return false;
    }
  }
  return true;
}

}  // namespace matching
}  // namespace openMVG
