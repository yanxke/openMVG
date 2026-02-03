#!/usr/bin/env python3
"""
Populate max translation distance constraints in sfm_data.json based on EXIF timestamps.

This script:
1. Reads an sfm_data.json file
2. Extracts EXIF timestamps (millisecond precision) from each image
3. For each pair of views that are temporally adjacent, computes:
   max_distance = max_velocity * delta_time
4. Adds these constraints to sfm_data.json as "max_translation_distance"

Usage:
    python populate_translation_constraints.py \
        --sfm_data sfm_data.json \
        --output sfm_data_with_constraints.json \
        --max_velocity 1.0 \
        --time_window 10.0

Arguments:
    --sfm_data: Input sfm_data.json file
    --output: Output sfm_data.json file with constraints
    --max_velocity: Maximum velocity in meters/second (default: 1.0)
    --time_window: Maximum time window in seconds to constrain (default: 10.0)
                   Pairs with dt > time_window will not get a constraint
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import exifread
except ImportError:
    print("ERROR: exifread not installed. Run: pip install exifread")
    exit(1)


def parse_exif_datetime(exif_tags) -> Optional[float]:
    """
    Extract timestamp from EXIF tags and return as seconds since epoch.
    
    Tries multiple EXIF fields in order:
    1. EXIF DateTimeOriginal + SubSecTimeOriginal
    2. EXIF DateTimeDigitized + SubSecTimeDigitized  
    3. Image DateTime + SubSecTime
    
    Returns:
        Timestamp in seconds since epoch (with millisecond precision), or None if not found
    """
    # Standard EXIF tag names for datetime and their corresponding sub-second tags
    datetime_pairs = [
        ('EXIF DateTimeOriginal', 'EXIF SubSecTimeOriginal'),
        ('EXIF DateTimeDigitized', 'EXIF SubSecTimeDigitized'),
        ('Image DateTime', 'EXIF SubSecTime'), # Sometimes it's here
        ('Image DateTime', 'Image SubSecTime'),
    ]
    
    for dt_tag, ss_tag in datetime_pairs:
        if dt_tag in exif_tags:
            dt_str = str(exif_tags[dt_tag])
            try:
                # Format: "YYYY:MM:DD HH:MM:SS"
                dt = datetime.strptime(dt_str, '%Y:%m:%d %H:%M:%S')
                
                # Default to 0 milliseconds
                ms = 0.0
                
                # Check for sub-seconds
                # Try the preferred pair first
                if ss_tag in exif_tags:
                    ss_val = str(exif_tags[ss_tag]).strip()
                    if ss_val:
                        ms = float("0." + ss_val)
                else:
                    # Fallback: Search for ANY tag containing 'SubSec'
                    for key in exif_tags.keys():
                        if 'SubSec' in key and ss_tag.split()[-1] in key:
                             ss_val = str(exif_tags[key]).strip()
                             if ss_val and ss_val.isdigit():
                                 ms = float("0." + ss_val)
                                 break
                
                return dt.timestamp() + ms
                
            except (ValueError, ZeroDivisionError):
                continue
    
    return None


def extract_timestamps_from_images(sfm_data: dict, root_path: str, image_dir: Optional[str] = None) -> Dict[int, float]:
    """
    Extract timestamps from images referenced in sfm_data.
    
    Args:
        sfm_data: Parsed sfm_data.json
        root_path: Root path for images (from sfm_data)
        image_dir: Manual override for image directory (host path)
        
    Returns:
        Dictionary mapping view_id -> timestamp (seconds since epoch)
    """
    timestamps = {}
    views = sfm_data.get('views', {})
    
    print(f"Extracting timestamps from {len(views)} views...")
    if image_dir:
        print(f"Image override directory: {image_dir}")
    
    for view_entry in views:
        if isinstance(view_entry, dict):
            view_id = view_entry['value']['ptr_wrapper']['data']['id_view']
            local_path = view_entry['value']['ptr_wrapper']['data'].get('local_path', '')
            filename = view_entry['value']['ptr_wrapper']['data'].get('filename', '')
        else:
            # Handle different JSON structures
            continue
            
        # Try to find the image
        img_path = None
        
        # 1. Try override directory if provided
        if image_dir:
            # Try directly in image_dir
            candidate = os.path.join(image_dir, filename)
            if os.path.exists(candidate):
                img_path = candidate
            else:
                # Try with the last component of local_path (e.g. if local_path is /data/orig_images)
                # but might be redundant if filename is unique.
                # Many SfM pipelines use the basename of the filename.
                candidate = os.path.join(image_dir, os.path.basename(filename))
                if os.path.exists(candidate):
                    img_path = candidate

        # 2. Try original path calculated from sfm_data + root_path
        if not img_path:
            if local_path:
                img_path = os.path.join(root_path, local_path, filename)
            else:
                img_path = os.path.join(root_path, filename)
            
        if not img_path or not os.path.exists(img_path):
            print(f"  Warning: Image not found: {img_path if img_path else filename}")
            continue
            
        # Read EXIF
        try:
            with open(img_path, 'rb') as f:
                tags = exifread.process_file(f, details=False)
                timestamp = parse_exif_datetime(tags)
                
                if timestamp is not None:
                    timestamps[view_id] = timestamp
                    dt = datetime.fromtimestamp(timestamp)
                    print(f"  View {view_id}: {dt.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]} ({os.path.basename(img_path)})")
                else:
                    print(f"  Warning: No timestamp found in EXIF for view {view_id}: {os.path.basename(img_path)}")
        except Exception as e:
            print(f"  Error reading EXIF from {img_path}: {e}")
            
    print(f"Successfully extracted {len(timestamps)} timestamps from {len(views)} views")
    return timestamps


def compute_distance_constraints(
    timestamps: Dict[int, float],
    max_velocity: float,
    time_window: float
) -> List[Tuple[int, int, float]]:
    """
    Compute pairwise distance constraints based on timestamps.
    
    Args:
        timestamps: Map of view_id -> timestamp
        max_velocity: Maximum velocity in m/s
        time_window: Maximum time window in seconds
        
    Returns:
        List of (view_i, view_j, max_distance) tuples
    """
    constraints = []
    
    # Sort views by timestamp
    sorted_views = sorted(timestamps.items(), key=lambda x: x[1])
    
    print(f"\nComputing distance constraints (max_velocity={max_velocity} m/s, time_window={time_window}s)...")
    
    for i in range(len(sorted_views)):
        view_i, time_i = sorted_views[i]
        
        # Look ahead in time to find nearby views
        for j in range(i + 1, len(sorted_views)):
            view_j, time_j = sorted_views[j]
            
            dt = abs(time_j - time_i)
            
            # Stop if time difference exceeds window
            if dt > time_window:
                break
                
            # Compute max distance
            max_distance = max_velocity * dt
            
            # Store constraint (always use lower ID first for consistency)
            constraints.append((min(view_i, view_j), max(view_i, view_j), max_distance))
    
    print(f"Generated {len(constraints)} distance constraints")
    
    # Print statistics
    if constraints:
        distances = [c[2] for c in constraints]
        print(f"  Distance range: {min(distances):.3f} - {max(distances):.3f} meters")
        print(f"  Mean distance: {sum(distances)/len(distances):.3f} meters")
    
    return constraints


def add_constraints_to_sfm_data(
    sfm_data: dict,
    constraints: List[Tuple[int, int, float]]
) -> dict:
    """
    Add max_translation_distance constraints to sfm_data.
    
    Args:
        sfm_data: Parsed sfm_data.json
        constraints: List of (view_i, view_j, max_distance) tuples
        
    Returns:
        Modified sfm_data dictionary
    """
    # Create the max_translation_distance field
    # Format: list of {"key": {"first": i, "second": j}, "value": distance}
    constraint_list = []
    
    for view_i, view_j, max_dist in constraints:
        constraint_list.append({
            "key": {
                "first": view_i,
                "second": view_j
            },
            "value": max_dist
        })
    
    sfm_data['max_translation_distance'] = constraint_list
    
    return sfm_data


def main():
    parser = argparse.ArgumentParser(
        description="Populate max translation distance constraints from EXIF timestamps"
    )
    parser.add_argument(
        '-p', '--project_dir',
        required=True,
        help='Project directory containing openmvg/matches/sfm_data.json and orig_images/'
    )
    parser.add_argument(
        '--max_velocity',
        type=float,
        default=1.0,
        help='Maximum velocity in meters/second (default: 1.0)'
    )
    parser.add_argument(
        '--time_window',
        type=float,
        default=10.0,
        help='Maximum time window in seconds for constraints (default: 10.0)'
    )
    
    args = parser.parse_args()

    # Define paths
    sfm_data_path = os.path.join(args.project_dir, 'openmvg/matches/sfm_data.json')
    sfm_data_orig_path = sfm_data_path + '.orig'
    image_dir = os.path.join(args.project_dir, 'orig_images')

    # Backup and File Choice Logic
    if not os.path.exists(sfm_data_path):
        print(f"ERROR: Base sfm_data.json not found at {sfm_data_path}")
        return 1

    if not os.path.exists(sfm_data_orig_path):
        import shutil
        print(f"Creating backup: {sfm_data_path} -> {sfm_data_orig_path}")
        shutil.copy2(sfm_data_path, sfm_data_orig_path)
        input_path = sfm_data_orig_path
    else:
        print(f"Backup already exists. Reading source from: {sfm_data_orig_path}")
        input_path = sfm_data_orig_path

    print(f"Reading SfM data: {input_path}")
    print(f"Writing result to: {sfm_data_path}")
    
    # Load sfm_data
    with open(input_path, 'r') as f:
        sfm_data = json.load(f)
    
    # Get root path from sfm_data (fallback for internal lookups)
    root_path = sfm_data.get('root_path', os.path.dirname(input_path))
    
    # Extract timestamps
    timestamps = extract_timestamps_from_images(sfm_data, root_path, image_dir)
    
    if not timestamps:
        print("ERROR: No timestamps found in EXIF data!")
        return 1
    
    # Compute constraints
    constraints = compute_distance_constraints(timestamps, args.max_velocity, args.time_window)
    
    if not constraints:
        print("WARNING: No constraints generated!")
        return 1
    
    # Add to sfm_data
    sfm_data = add_constraints_to_sfm_data(sfm_data, constraints)
    
    # Save output (always overwrite the main file)
    with open(sfm_data_path, 'w') as f:
        json.dump(sfm_data, f, indent=2)
    
    print(f"Successfully updated {sfm_data_path} with translation constraints.")
    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
