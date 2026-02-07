#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/motion_from_essential.hpp"
#include "openMVG/multiview/solver_essential_eight_point.hpp"
#include "openMVG/numeric/numeric.h"
#include "openMVG/sfm/pipelines/sfm_features_provider.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/matching_image_collection/GeometricFilter.hpp"
#include "openMVG/matching_image_collection/E_ACRobust.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <cmath>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::sfm;
using namespace openMVG::matching_image_collection;

int main(int argc, char **argv) {
    CmdLine cmd;

    std::string sSfM_Data_Filename;
    std::string sMatchesFilename;
    std::string sOutputMatchesFilename;
    std::string sHeadingsFilename;
    std::string sFeatDir;

    double yaw_tol = 15.0;   // degrees
    double pitch_tol = 10.0; // degrees
    double roll_tol = 10.0;  // degrees
    double alt_tol = 0.2;    // max |ty| in unit t (approx sin(angle))

    cmd.add( make_option('i', sSfM_Data_Filename, "input_file") );
    cmd.add( make_option('m', sMatchesFilename, "matches") );
    cmd.add( make_option('o', sOutputMatchesFilename, "output_file") );
    cmd.add( make_option('p', sHeadingsFilename, "headings") );
    cmd.add( make_option('f', sFeatDir, "feat_dir") );
    cmd.add( make_option('y', yaw_tol, "yaw_tol") );
    cmd.add( make_option('x', pitch_tol, "pitch_tol") );
    cmd.add( make_option('z', roll_tol, "roll_tol") );
    cmd.add( make_option('a', alt_tol, "alt_tol") );

    try {
        if (argc == 1) throw std::string("Invalid parameter.");
        cmd.process(argc, argv);
    } catch (const std::string& s) {
        std::cerr << "Usage: " << argv[0] << " -i sfm_data.json -m matches.e.bin -o matches.e.filtered.bin -p headings.txt -f feat_dir\n";
        return EXIT_FAILURE;
    }

    // 1. Load SfM Data
    SfM_Data sfm_data;
    if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(VIEWS|INTRINSICS))) {
        std::cerr << "Cannot load sfm_data\n";
        return EXIT_FAILURE;
    }

    // 2. Load Matches
    PairWiseMatches map_Matches;
    if (!Load(map_Matches, sMatchesFilename)) {
        std::cerr << "Cannot load matches\n";
        return EXIT_FAILURE;
    }

    // 3. (Skipped) Headings are read directly from EXIF

    // 4. Features Provider
    std::shared_ptr<Features_Provider> feats_provider = std::make_shared<Features_Provider>();
    const std::string sImageDescriber = stlplus::create_filespec(sFeatDir, "image_describer", "json");
    std::unique_ptr<openMVG::features::Regions> regions_type = openMVG::features::Init_region_type_from_file(sImageDescriber);
    if (!regions_type || !feats_provider->load(sfm_data, sFeatDir, regions_type)) {
        std::cerr << "Cannot load features from: " << sFeatDir << "\n";
        return EXIT_FAILURE;
    }

    // Cache headings
    std::map<IndexT, double> map_headings;
    for (const auto & view_ptr : sfm_data.GetViews()) {
        std::unique_ptr<openMVG::exif::Exif_IO> exifIO(new openMVG::exif::Exif_IO_EasyExif(stlplus::folder_append_separator(sfm_data.s_root_path) + view_ptr.second->s_Img_path));
        double heading;
        if (exifIO->GPSImageDirection(&heading)) {
            map_headings[view_ptr.first] = heading;
        }
    }
    std::cout << "Loaded headings for " << map_headings.size() << " images.\n";

    PairWiseMatches filtered_Matches;
    int count_total = map_Matches.size();
    int count_kept = 0;
    int count_fail_pitch = 0;
    int count_fail_roll = 0;
    int count_fail_yaw = 0;
    int count_fail_alt = 0;

    std::vector<double> hist_pitch, hist_roll, hist_yaw_err, hist_alt;

    for (const auto & pair_it : map_Matches) {
        IndexT i = pair_it.first.first;
        IndexT j = pair_it.first.second;
        const IndMatches & matches = pair_it.second;

        if (matches.size() < 8) continue; 

        // Get intrinsics
        auto it_i = sfm_data.GetIntrinsics().find(sfm_data.GetViews().at(i)->id_intrinsic);
        auto it_j = sfm_data.GetIntrinsics().find(sfm_data.GetViews().at(j)->id_intrinsic);
        if (it_i == sfm_data.GetIntrinsics().end() || it_j == sfm_data.GetIntrinsics().end()) continue;

        const cameras::IntrinsicBase * intrinsic1 = it_i->second.get();
        const cameras::IntrinsicBase * intrinsic2 = it_j->second.get();

        const auto & feat1 = feats_provider->feats_per_view.at(i);
        const auto & feat2 = feats_provider->feats_per_view.at(j);

        Mat3X bearing_x1(3, matches.size());
        Mat3X bearing_x2(3, matches.size());
        
        int col_idx = 0;
        for (const auto & m : matches) {
            Vec2 p1 = feat1[m.i_].coords().template cast<double>();
            Vec2 p2 = feat2[m.j_].coords().template cast<double>();
            
            Vec3 b1 = (intrinsic1->ima2cam(p1)).homogeneous().normalized();
            Vec3 b2 = (intrinsic2->ima2cam(p2)).homogeneous().normalized();
            
            bearing_x1.col(col_idx) = b1;
            bearing_x2.col(col_idx) = b2;
            col_idx++;
        }

        // We assume matches are already geometrically filtered (inliers).
        // 5-point algorithm (or 8-point) to recover E, then R, t.
        if (bearing_x1.cols() >= 8) {
            std::vector<Mat3> Es;
            EightPointRelativePoseSolver::Solve(bearing_x1, bearing_x2, &Es);
            if(Es.empty()) continue;
            Mat3 E = Es[0];
            
            geometry::Pose3 pose;
            std::vector<uint32_t> inliers_indices(bearing_x1.cols());
            std::iota(inliers_indices.begin(), inliers_indices.end(), 0);
            
            std::vector<uint32_t> selected_points;
            if (RelativePoseFromEssential(bearing_x1, bearing_x2, E, inliers_indices, &pose, &selected_points)) {
                Mat3 R = pose.rotation();
                Vec3 t = pose.translation();

                // Decompose R to Pitch, Yaw, Roll
                // Camera coords: x-right, y-down, z-forward
                // Pitch (around x), Yaw (around y), Roll (around z)
                double pitch = asin(clamp(R(1, 2), -1.0, 1.0));
                double yaw = atan2(-R(0, 2), R(2, 2));
                double roll = atan2(-R(1, 0), R(1, 1));

                pitch *= 180.0 / M_PI;
                yaw *= 180.0 / M_PI;
                roll *= 180.0 / M_PI;

                bool fail_p = std::abs(pitch) > pitch_tol;
                bool fail_r = std::abs(roll) > roll_tol;
                bool fail_a = std::abs(t(1)) > alt_tol;
                bool fail_y = false;
                double diff = 0.0;

                if (map_headings.count(i) && map_headings.count(j)) {
                    double h_i = map_headings[i];
                    double h_j = map_headings[j];
                    double dyaw_expected = h_j - h_i;
                    while (dyaw_expected > 180) dyaw_expected -= 360;
                    while (dyaw_expected < -180) dyaw_expected += 360;
                    
                    diff = std::abs(yaw - dyaw_expected);
                    if (diff > 180) diff = 360 - diff;
                    if (diff > yaw_tol) fail_y = true;
                }

                if (fail_p) count_fail_pitch++;
                if (fail_r) count_fail_roll++;
                if (fail_a) count_fail_alt++;
                if (fail_y) count_fail_yaw++;

                if (!fail_p && !fail_r && !fail_a && !fail_y) {
                    filtered_Matches[pair_it.first] = matches;
                    count_kept++;
                    hist_pitch.push_back(std::abs(pitch));
                    hist_roll.push_back(std::abs(roll));
                    hist_alt.push_back(std::abs(t(1)));
                    if (map_headings.count(i) && map_headings.count(j)) {
                        hist_yaw_err.push_back(diff);
                    }
                }
            }
        }
    }

    std::cout << "FilterByPriors: Kept " << count_kept << " / " << count_total << " pairs.\n";
    std::cout << "  Fail Pitch: " << count_fail_pitch << "\n";
    std::cout << "  Fail Roll:  " << count_fail_roll << "\n";
    std::cout << "  Fail Alt:   " << count_fail_alt << "\n";
    std::cout << "  Fail Yaw:   " << count_fail_yaw << "\n";

    auto print_hist = [](const std::vector<double>& v, const std::string& name, double max_val) {
        if (v.empty()) return;
        int n_bins = 10;
        std::vector<int> bins(n_bins, 0);
        for (double val : v) {
            int idx = std::min(n_bins - 1, (int)(val / max_val * n_bins));
            if (idx >= 0) bins[idx]++;
        }
        std::cout << "\nHistogram for " << name << " (0 to " << max_val << "):\n";
        double step = max_val / n_bins;
        for (int i = 0; i < n_bins; ++i) {
            std::cout << "[" << i*step << "-" << (i+1)*step << "]: " << bins[i] << "\n";
        }
    };

    print_hist(hist_pitch, "Pitch Error (deg)", pitch_tol);
    print_hist(hist_roll, "Roll Error (deg)", roll_tol);
    print_hist(hist_alt, "Altitude Error (unit t)", alt_tol);
    print_hist(hist_yaw_err, "Yaw Mismatch (deg)", yaw_tol);

    if (!Save(filtered_Matches, sOutputMatchesFilename)) {
        std::cerr << "Cannot save filtered matches\n";
        return EXIT_FAILURE;
    }

    // Also save as JSON for easier parsing
    const std::string sJsonFilename = stlplus::create_filespec(
      stlplus::folder_part(sOutputMatchesFilename),
      stlplus::basename_part(sOutputMatchesFilename),
      "json");
    if ( !SaveJson( filtered_Matches, sJsonFilename ) )
    {
      std::cerr << "Warning: Cannot save JSON matches to: " << sJsonFilename << "\n";
    }
    else
    {
      std::cout << "Saved JSON matches to: " << sJsonFilename << "\n";
    }

    return EXIT_SUCCESS;
}
