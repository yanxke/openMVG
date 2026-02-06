// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/cameras.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"
#include "openMVG/exif/sensor_width_database/ParseDatabase.hpp"
#include "openMVG/geodesy/geodesy.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"
#include "openMVG/sfm/sfm_data_utils.hpp"
#include "openMVG/sfm/sfm_view.hpp"
#include "openMVG/sfm/sfm_view_priors.hpp"
#include "openMVG/system/loggerprogress.hpp"
#include "openMVG/types.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace openMVG;
using namespace openMVG::cameras;
using namespace openMVG::exif;
using namespace openMVG::geodesy;
using namespace openMVG::image;
using namespace openMVG::sfm;

/// Parse IMU rotation matrix from UserComment EXIF field
/// Format: "Rotation: r11,r12,r13,r21,r22,r23,r31,r32,r33"
/// Returns Device->World rotation in ENU frame
bool ParseImuRotationFromUserComment(const std::string & comment, openMVG::Mat3 & rotation_dw)
{
  const std::string key = "Rotation:";
  const std::size_t pos = comment.find(key);
  if (pos == std::string::npos)
    return false;

  std::string rot_part = comment.substr(pos + key.size());
  for (char & c : rot_part)
  {
    if (c == '\n' || c == '\r')
      c = ' ';
  }

  std::vector<double> vals;
  vals.reserve(9);
  std::string token;
  std::stringstream ss(rot_part);
  while (std::getline(ss, token, ','))
  {
    std::stringstream t(token);
    double v = 0.0;
    if (t >> v)
      vals.push_back(v);
  }
  if (vals.size() < 9)
    return false;

  // Validate all values are finite (not NaN or Inf)
  for (const double& val : vals)
  {
    if (!std::isfinite(val))
    {
      std::cerr << "Warning: IMU rotation contains non-finite value, ignoring" << std::endl;
      return false;
    }
  }

  rotation_dw << vals[0], vals[1], vals[2],
                  vals[3], vals[4], vals[5],
                  vals[6], vals[7], vals[8];
  return true;
}

/// Parse XMP stepsSinceTaskStart from image file
/// XMP namespace: https://clobotics.com/storecapture/1.0/
bool ParseXmpStepCounter(const std::string & filename, int & step_counter)
{
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open())
    return false;

  // Read file content (limit to first 1MB for XMP search)
  const size_t max_read = 1024 * 1024;
  std::vector<char> buffer(max_read);
  file.read(buffer.data(), max_read);
  std::streamsize bytes_read = file.gcount();
  file.close();

  if (bytes_read <= 0)
    return false;

  std::string content(buffer.data(), bytes_read);

  // Look for stepsSinceTaskStart in XMP data
  // Format: sc:stepsSinceTaskStart="1062"
  const std::string key = "stepsSinceTaskStart";
  size_t pos = content.find(key);

  if (pos == std::string::npos)
    return false;

  // Find the value after the key (format: stepsSinceTaskStart="VALUE")
  // Move to the end of the key
  pos += key.length();

  // Find the opening quote
  pos = content.find('"', pos);
  if (pos == std::string::npos)
    return false;
  pos++; // Move past opening quote

  // Find the closing quote
  size_t end_pos = content.find('"', pos);
  if (end_pos == std::string::npos)
    return false;

  std::string value_str = content.substr(pos, end_pos - pos);

  // Trim whitespace
  value_str.erase(0, value_str.find_first_not_of(" \t\n\r"));
  value_str.erase(value_str.find_last_not_of(" \t\n\r") + 1);

  // Parse integer
  try
  {
    step_counter = std::stoi(value_str);
    return true;
  }
  catch (...)
  {
    return false;
  }
}

/// Check that Kmatrix is a string like "f;0;ppx;0;f;ppy;0;0;1"
/// With f,ppx,ppy as valid numerical value
bool checkIntrinsicStringValidity(const std::string & Kmatrix, double & focal, double & ppx, double & ppy)
{
  std::vector<std::string> vec_str;
  stl::split(Kmatrix, ';', vec_str);
  if (vec_str.size() != 9)  {
    OPENMVG_LOG_ERROR << "\n Missing ';' character";
    return false;
  }
  // Check that all K matrix value are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i) {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      OPENMVG_LOG_ERROR << "\n Used an invalid not a number character";
      return false;
    }
    if (i==0) focal = readvalue;
    if (i==2) ppx = readvalue;
    if (i==5) ppy = readvalue;
  }
  return true;
}

bool getGPS
(
  const std::string & filename,
  const int & GPS_to_XYZ_method,
  Vec3 & pose_center
)
{
  std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
  if (exifReader)
  {
    // Try to parse EXIF metada & check existence of EXIF data
    if ( exifReader->open( filename ) && exifReader->doesHaveExifInfo() )
    {
      // Check existence of GPS coordinates
      double latitude, longitude, altitude;
      if ( exifReader->GPSLatitude( &latitude ) &&
           exifReader->GPSLongitude( &longitude ) &&
           exifReader->GPSAltitude( &altitude ) )
      {
        // Add ECEF or UTM XYZ position to the GPS position array
        switch (GPS_to_XYZ_method)
        {
          case 1:
            pose_center = lla_to_utm( latitude, longitude, altitude );
            break;
          case 0:
          default:
            pose_center = lla_to_ecef( latitude, longitude, altitude );
            break;
        }
        return true;
      }
    }
  }
  return false;
}

bool computeSensorWidthFromExif(const Exif_IO &exifReader,
  const double image_width_px,
  const double image_height_px,
  double &sensor_width_mm,
  std::string &method)
{
  const double focal_mm = static_cast<double>(exifReader.getFocal());
  const double focal_35mm = static_cast<double>(exifReader.getFocalLengthIn35mm());
  if (focal_mm <= 0.0 || focal_35mm <= 0.0)
  {
    // Fall through to focal plane resolution-based estimate.
  }
  else
  {
    sensor_width_mm = 36.0 * focal_mm / focal_35mm;
    method = "35mm_equivalent";
    return sensor_width_mm > 0.0;
  }

  const double x_resolution = static_cast<double>(exifReader.getFocalPlaneXResolution());
  const double y_resolution = static_cast<double>(exifReader.getFocalPlaneYResolution());
  const int resolution_unit = exifReader.getFocalPlaneResolutionUnit();
  double unit_mm = 0.0;
  switch (resolution_unit)
  {
    case 2: // inch
      unit_mm = 25.4;
      break;
    case 3: // centimeter
      unit_mm = 10.0;
      break;
    case 4: // millimeter
      unit_mm = 1.0;
      break;
    case 5: // micrometer
      unit_mm = 0.001;
      break;
    default:
      unit_mm = 0.0;
      break;
  }

  if (unit_mm <= 0.0)
  {
    return false;
  }

  if (x_resolution > 0.0 && image_width_px > 0.0)
  {
    sensor_width_mm = image_width_px / x_resolution * unit_mm;
    method = "focal_plane_x_resolution";
    return sensor_width_mm > 0.0;
  }

  if (y_resolution > 0.0 && image_height_px > 0.0)
  {
    sensor_width_mm = image_height_px / y_resolution * unit_mm;
    method = "focal_plane_y_resolution";
    return sensor_width_mm > 0.0;
  }

  return false;
}

/// Check string of prior weights
std::pair<bool, Vec3> checkPriorWeightsString
(
  const std::string &sWeights
)
{
  std::pair<bool, Vec3> val(true, Vec3::Zero());
  std::vector<std::string> vec_str;
  stl::split(sWeights, ';', vec_str);
  if (vec_str.size() != 3)
  {
    OPENMVG_LOG_ERROR << "Missing ';' character in prior weights";
    val.first = false;
  }
  // Check that all weight values are valid numbers
  for (size_t i = 0; i < vec_str.size(); ++i)
  {
    double readvalue = 0.0;
    std::stringstream ss;
    ss.str(vec_str[i]);
    if (! (ss >> readvalue) )  {
      OPENMVG_LOG_ERROR << "Used an invalid not a number character in local frame origin";
      val.first = false;
    }
    val.second[i] = readvalue;
  }
  return val;
}
//
// Create the description of an input image dataset for OpenMVG toolsuite
// - Export a SfM_Data file with View & Intrinsic data
//
int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImageDir,
    sfileDatabase = "",
    sOutputDir = "",
    sKmatrix;

  std::string sPriorWeights = "1.0;1.0;1.0";
  std::pair<bool, Vec3> prior_w_info(false, Vec3());

  int i_User_camera_model = PINHOLE_CAMERA_RADIAL3;

  bool b_Group_camera_model = true;

  int i_GPS_XYZ_method = 0;

  double focal_pixels = -1.0;

  cmd.add( make_option('i', sImageDir, "imageDirectory") );
  cmd.add( make_option('d', sfileDatabase, "sensorWidthDatabase") );
  cmd.add( make_option('o', sOutputDir, "outputDirectory") );
  cmd.add( make_option('f', focal_pixels, "focal") );
  cmd.add( make_option('k', sKmatrix, "intrinsics") );
  cmd.add( make_option('c', i_User_camera_model, "camera_model") );
  cmd.add( make_option('g', b_Group_camera_model, "group_camera_model") );
  cmd.add( make_switch('P', "use_pose_prior") );
  cmd.add( make_option('W', sPriorWeights, "prior_weights"));
  cmd.add( make_option('m', i_GPS_XYZ_method, "gps_to_xyz_method") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch (const std::string& s) {
    OPENMVG_LOG_INFO << "Usage: " << argv[0] << '\n'
      << "[-i|--imageDirectory]\n"
      << "[-d|--sensorWidthDatabase]\n"
      << "[-o|--outputDirectory]\n"
      << "[-f|--focal] (pixels)\n"
      << "[-k|--intrinsics] Kmatrix: \"f;0;ppx;0;f;ppy;0;0;1\"\n"
      << "[-c|--camera_model] Camera model type:\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA) << ": Pinhole\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_RADIAL1) << ": Pinhole radial 1\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_RADIAL3) << ": Pinhole radial 3 (default)\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_BROWN) << ": Pinhole brown 2\n"
      << "\t" << static_cast<int>(PINHOLE_CAMERA_FISHEYE) << ": Pinhole with a simple Fish-eye distortion\n"
      << "\t" << static_cast<int>(CAMERA_SPHERICAL) << ": Spherical camera\n"
      << "[-g|--group_camera_model]\n"
      << "\t 0-> each view have it's own camera intrinsic parameters,\n"
      << "\t 1-> (default) view can share some camera intrinsic parameters\n"
      << "\n"
      << "[-P|--use_pose_prior] Use pose prior if GPS EXIF pose is available"
      << "[-W|--prior_weights] \"x;y;z;\" of weights for each dimension of the prior (default: 1.0)\n"
      << "[-m|--gps_to_xyz_method] XZY Coordinate system:\n"
      << "\t 0: ECEF (default)\n"
      << "\t 1: UTM";

      OPENMVG_LOG_ERROR << s;
      return EXIT_FAILURE;
  }

  const bool b_Use_pose_prior = cmd.used('P');
  OPENMVG_LOG_INFO << " You called : " << argv[0]
    << "\n--imageDirectory " << sImageDir
    << "\n--sensorWidthDatabase " << sfileDatabase
    << "\n--outputDirectory " << sOutputDir
    << "\n--focal " << focal_pixels
    << "\n--intrinsics " << sKmatrix
    << "\n--camera_model " << i_User_camera_model
    << "\n--group_camera_model " << b_Group_camera_model
    << "\n--use_pose_prior " << b_Use_pose_prior
    << "\n--prior_weights " << sPriorWeights
    << "\n--gps_to_xyz_method " << i_GPS_XYZ_method;

  // Expected properties for each image
  double width = -1, height = -1, focal = -1, ppx = -1,  ppy = -1;

  const EINTRINSIC e_User_camera_model = EINTRINSIC(i_User_camera_model);

  if ( !stlplus::folder_exists( sImageDir ) )
  {
    OPENMVG_LOG_ERROR << "The input directory doesn't exist";
    return EXIT_FAILURE;
  }

  if (sOutputDir.empty())
  {
    OPENMVG_LOG_ERROR << "Invalid output directory";
    return EXIT_FAILURE;
  }

  if ( !stlplus::folder_exists( sOutputDir ) )
  {
    if ( !stlplus::folder_create( sOutputDir ))
    {
      OPENMVG_LOG_ERROR << "Cannot create output directory";
      return EXIT_FAILURE;
    }
  }

  const std::string sOutputSfMData =
    stlplus::create_filespec( sOutputDir, "sfm_data.json" );
  if ( stlplus::file_exists( sOutputSfMData ) )
  {
    OPENMVG_LOG_INFO
      << "sfm_data.json already exists, skipping image listing: "
      << sOutputSfMData;
    return EXIT_SUCCESS;
  }

  if (sKmatrix.size() > 0 &&
    !checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy) )
  {
    OPENMVG_LOG_ERROR << "Invalid K matrix input";
    return EXIT_FAILURE;
  }

  if (sKmatrix.size() > 0 && focal_pixels != -1.0)
  {
    OPENMVG_LOG_ERROR << "Cannot combine -f and -k options";
    return EXIT_FAILURE;
  }

  std::vector<Datasheet> vec_database;
  if (!sfileDatabase.empty())
  {
    if ( !parseDatabase( sfileDatabase, vec_database ) )
    {
      OPENMVG_LOG_ERROR
       << "Invalid input database: " << sfileDatabase
       << ", please specify a valid file.";
      return EXIT_FAILURE;
    }
  }

  // Check if prior weights are given
  if (b_Use_pose_prior)
  {
    prior_w_info = checkPriorWeightsString(sPriorWeights);
  }

  std::vector<std::string> vec_image = stlplus::folder_files( sImageDir );
  std::sort(vec_image.begin(), vec_image.end());

  // Configure an empty scene with Views and their corresponding cameras
  SfM_Data sfm_data;
  sfm_data.s_root_path = sImageDir; // Setup main image root_path
  Views & views = sfm_data.views;
  Intrinsics & intrinsics = sfm_data.intrinsics;

  bool sensor_report_printed = false;

  system::LoggerProgress my_progress_bar(vec_image.size(), "- Listing images -" );
  std::ostringstream error_report_stream;
  for ( std::vector<std::string>::const_iterator iter_image = vec_image.begin();
    iter_image != vec_image.end();
    ++iter_image, ++my_progress_bar )
  {
    // Read meta data to fill camera parameter (w,h,focal,ppx,ppy) fields.
    width = height = ppx = ppy = focal = -1.0;

    const std::string sImageFilename = stlplus::create_filespec( sImageDir, *iter_image );
    const std::string sImFilenamePart = stlplus::filename_part(sImageFilename);

    // Test if the image format is supported:
    if (openMVG::image::GetFormat(sImageFilename.c_str()) == openMVG::image::Unknown)
    {
      error_report_stream
          << sImFilenamePart << ": Unkown image file format." << "\n";
      continue; // image cannot be opened
    }

    if (sImFilenamePart.find("mask.png") != std::string::npos
       || sImFilenamePart.find("_mask.png") != std::string::npos)
    {
      error_report_stream
          << sImFilenamePart << " is a mask image" << "\n";
      continue;
    }

    ImageHeader imgHeader;
    if (!openMVG::image::ReadImageHeader(sImageFilename.c_str(), &imgHeader))
      continue; // image cannot be read

    width = imgHeader.width;
    height = imgHeader.height;
    ppx = width / 2.0;
    ppy = height / 2.0;


    // Consider the case where the focal is provided manually
    if (sKmatrix.size() > 0) // Known user calibration K matrix
    {
      if (!checkIntrinsicStringValidity(sKmatrix, focal, ppx, ppy))
        focal = -1.0;
    }
    else // User provided focal length value
      if (focal_pixels != -1 )
        focal = focal_pixels;

    // If not manually provided or wrongly provided
    if (focal == -1)
    {
      std::unique_ptr<Exif_IO> exifReader(new Exif_IO_EasyExif);
      exifReader->open( sImageFilename );

      const bool bHaveValidExifMetadata =
        exifReader->doesHaveExifInfo()
        && !exifReader->getModel().empty()
        && !exifReader->getBrand().empty();

      if (bHaveValidExifMetadata) // If image contains meta data
      {
        // Handle case where focal length is equal to 0
        if (exifReader->getFocal() == 0.0f)
        {
          error_report_stream
            << stlplus::basename_part(sImageFilename) << ": Focal length is missing." << "\n";
          focal = -1.0;
        }
        else
        // Create the image entry in the list file
        {
          const std::string sCamModel = exifReader->getBrand() + " " + exifReader->getModel();

          Datasheet datasheet;
          const bool has_db_entry = getInfo( sCamModel, vec_database, datasheet );
          double exif_sensor_width_mm = -1.0;
          std::string exif_sensor_method;
          const bool has_exif_sensor_width = computeSensorWidthFromExif(
            *exifReader, width, height, exif_sensor_width_mm, exif_sensor_method);

          if (!sensor_report_printed)
          {
            std::ostringstream report;
            report << "Sensor width report for model: " << sCamModel;
            if (has_exif_sensor_width)
            {
            report << "\n- EXIF-derived sensor width (35mm equiv): " << exif_sensor_width_mm << " mm";
            if (!exif_sensor_method.empty())
            {
              report << " [" << exif_sensor_method << "]";
            }
          }
          else
          {
            report << "\n- EXIF-derived sensor width (35mm equiv): unavailable";
          }
            if (!vec_database.empty() && has_db_entry)
            {
              report << "\n- Sensor database width: " << datasheet.sensorSize_ << " mm";
            }
            else if (!vec_database.empty())
            {
              report << "\n- Sensor database width: not found for model";
            }
            OPENMVG_LOG_INFO << report.str();
            sensor_report_printed = true;
          }

          if (has_db_entry)
          {
            // The camera model was found in the database so we can compute it's approximated focal length
            const double ccdw = datasheet.sensorSize_;
            focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          }
          else if (has_exif_sensor_width)
          {
            const double ccdw = exif_sensor_width_mm;
            focal = std::max ( width, height ) * exifReader->getFocal() / ccdw;
          }
          else
          {
            error_report_stream
              << stlplus::basename_part(sImageFilename)
              << "\" model \"" << sCamModel << "\" doesn't exist in the database" << "\n"
              << "Please consider add your camera model and sensor width in the database." << "\n";
          }
        }
      }
    }
    // Build intrinsic parameter related to the view
    std::shared_ptr<IntrinsicBase> intrinsic;

    if (focal > 0 && ppx > 0 && ppy > 0 && width > 0 && height > 0)
    {
      // Create the desired camera type
      switch (e_User_camera_model)
      {
        case PINHOLE_CAMERA:
          intrinsic = std::make_shared<Pinhole_Intrinsic>
            (width, height, focal, ppx, ppy);
        break;
        case PINHOLE_CAMERA_RADIAL1:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K1>
            (width, height, focal, ppx, ppy, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_RADIAL3:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Radial_K3>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0);  // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_BROWN:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Brown_T2>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case PINHOLE_CAMERA_FISHEYE:
          intrinsic = std::make_shared<Pinhole_Intrinsic_Fisheye>
            (width, height, focal, ppx, ppy, 0.0, 0.0, 0.0, 0.0); // setup no distortion as initial guess
        break;
        case CAMERA_SPHERICAL:
           intrinsic = std::make_shared<Intrinsic_Spherical>
             (width, height);
        break;
        default:
          OPENMVG_LOG_ERROR << "Error: unknown camera model: " << (int) e_User_camera_model;
          return EXIT_FAILURE;
      }
    }

    // Build the view corresponding to the image
    Vec3 pose_center;
    if (getGPS(sImageFilename, i_GPS_XYZ_method, pose_center) && b_Use_pose_prior)
    {
      ViewPriors v(*iter_image, views.size(), views.size(), views.size(), width, height);

      // Add intrinsic related to the image (if any)
      if (!intrinsic)
      {
        //Since the view have invalid intrinsic data
        // (export the view, with an invalid intrinsic field value)
        v.id_intrinsic = UndefinedIndexT;
      }
      else
      {
        // Add the defined intrinsic to the sfm_container
        intrinsics[v.id_intrinsic] = intrinsic;
      }

      v.b_use_pose_center_ = true;
      v.pose_center_ = pose_center;
      // prior weights
      if (prior_w_info.first == true)
      {
        v.center_weight_ = prior_w_info.second;
      }

      // Try to read GPS compass heading and IMU rotation
      std::unique_ptr<Exif_IO> exifReader_heading(new Exif_IO_EasyExif);
      if (exifReader_heading->open(sImageFilename))
      {
        double heading;
        if (exifReader_heading->GPSImageDirection(&heading))
        {
          v.b_has_heading_ = true;
          v.gps_heading_ = heading;
        }

        // Try to read IMU rotation from UserComment
        std::string user_comment;
        if (exifReader_heading->UserComment(&user_comment))
        {
          Mat3 R_dw;
          if (ParseImuRotationFromUserComment(user_comment, R_dw))
          {
            v.b_has_imu_rotation_ = true;
            v.imu_rotation_ = R_dw;
          }
        }
      }

      // Try to read XMP step counter
      int step_counter;
      if (ParseXmpStepCounter(sImageFilename, step_counter))
      {
        v.b_has_step_counter_ = true;
        v.step_counter_ = step_counter;
      }

      // Add the view to the sfm_container
      views[v.id_view] = std::make_shared<ViewPriors>(v);
    }
    else
    {
      // Check if we have GPS heading, IMU rotation, or XMP step counter even without GPS position
      std::unique_ptr<Exif_IO> exifReader_heading(new Exif_IO_EasyExif);
      double heading;
      bool has_heading = false;
      Mat3 R_dw;
      bool has_imu_rotation = false;
      int step_counter;
      bool has_step_counter = false;

      if (exifReader_heading->open(sImageFilename))
      {
        if (exifReader_heading->GPSImageDirection(&heading))
        {
          has_heading = true;
        }

        // Try to read IMU rotation from UserComment
        std::string user_comment;
        if (exifReader_heading->UserComment(&user_comment))
        {
          if (ParseImuRotationFromUserComment(user_comment, R_dw))
          {
            has_imu_rotation = true;
          }
        }
      }

      // Try to read XMP step counter
      if (ParseXmpStepCounter(sImageFilename, step_counter))
      {
        has_step_counter = true;
      }

      // If we have any metadata (heading, IMU, or step counter), use ViewPriors to store it
      if (has_heading || has_imu_rotation || has_step_counter)
      {
        ViewPriors v(*iter_image, views.size(), views.size(), views.size(), width, height);

        // Add intrinsic related to the image (if any)
        if (!intrinsic)
        {
          v.id_intrinsic = UndefinedIndexT;
        }
        else
        {
          intrinsics[v.id_intrinsic] = intrinsic;
        }

        if (has_heading)
        {
          v.b_has_heading_ = true;
          v.gps_heading_ = heading;
        }

        if (has_imu_rotation)
        {
          v.b_has_imu_rotation_ = true;
          v.imu_rotation_ = R_dw;
        }

        if (has_step_counter)
        {
          v.b_has_step_counter_ = true;
          v.step_counter_ = step_counter;
        }

        // Add the view to the sfm_container
        views[v.id_view] = std::make_shared<ViewPriors>(v);
      }
      else
      {
        // No heading available, use regular View
        View v(*iter_image, views.size(), views.size(), views.size(), width, height);

        // Add intrinsic related to the image (if any)
        if (!intrinsic)
        {
          //Since the view have invalid intrinsic data
          // (export the view, with an invalid intrinsic field value)
          v.id_intrinsic = UndefinedIndexT;
        }
        else
        {
          // Add the defined intrinsic to the sfm_container
          intrinsics[v.id_intrinsic] = intrinsic;
        }

        // Add the view to the sfm_container
        views[v.id_view] = std::make_shared<View>(v);
      }
    }
  }

  // Display saved warning & error messages if any.
  if (!error_report_stream.str().empty())
  {
    OPENMVG_LOG_WARNING
      << "Warning & Error messages:\n"
      << error_report_stream.str();
  }

  // Group camera that share common properties if desired (leads to more faster & stable BA).
  if (b_Group_camera_model)
  {
    GroupSharedIntrinsics(sfm_data);
  }

  // Store SfM_Data views & intrinsic data
  if (!Save(
    sfm_data,
    sOutputSfMData.c_str(),
    ESfM_Data(VIEWS|INTRINSICS)))
  {
    return EXIT_FAILURE;
  }

  OPENMVG_LOG_INFO
    << "SfMInit_ImageListing report:\n"
    << "listed #File(s): " << vec_image.size() << "\n"
    << "usable #File(s) listed in sfm_data: " << sfm_data.GetViews().size() << "\n"
    << "usable #Intrinsic(s) listed in sfm_data: " << sfm_data.GetIntrinsics().size();

  return EXIT_SUCCESS;
}
