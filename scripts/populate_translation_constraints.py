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
    # Try different EXIF fields in order of preference
    datetime_fields = [
        ('EXIF DateTimeOriginal', 'EXIF SubSecTimeOriginal'),
        ('EXIF DateTimeDigitized', 'EXIF SubSecTimeDigitized'),
        ('Image DateTime', 'Image SubSecTime'),
    ]
    
    for dt_field, subsec_field in datetime_fields:
        if dt_field in exif_tags:
            # Parse main datetime
            dt_str = str(exif_tags[dt_field])
            try:
                # Format: "YYYY:MM:DD HH:MM:SS"
                dt = datetime.strptime(dt_str, '%Y:%m:%d %H:%M:%S')
                
                # Add milliseconds if available
                milliseconds = 0
                if subsec_field in exif_tags:
                    subsec_str = str(exif_tags[subsec_field])
                    try:
                        # SubSec is usually a fraction like "123" meaning 0.123 seconds
                        # Pad or truncate to 3 digits for milliseconds
                        subsec_str = subsec_str.ljust(3, '0')[:3]
                        milliseconds = int(subsec_str)
                    except ValueError:
                        pass
                
                # Convert to timestamp
                timestamp = dt.timestamp() + milliseconds / 1000.0
                return timestamp
                
            except ValueError as e:
                print(f"  Warning: Could not parse datetime '{dt_str}': {e}")
                continue
    
    return None


def extract_timestamps_from_images(sfm_data: dict, root_path: str) -> Dict[int, float]:
    """
    Extract timestamps from images referenced in sfm_data.
    
    Args:
        sfm_data: Parsed sfm_data.json
        root_path: Root path for images
        
    Returns:
        Dictionary mapping view_id -> timestamp (seconds since epoch)
    """
    timestamps = {}
    views = sfm_data.get('views', {})
    
    print(f"Extracting timestamps from {len(views)} views...")
    
    for view_entry in views:
        if isinstance(view_entry, dict):
            view_id = view_entry['value']['ptr_wrapper']['data']['id_view']
            local_path = view_entry['value']['ptr_wrapper']['data'].get('local_path', '')
            filename = view_entry['value']['ptr_wrapper']['data'].get('filename', '')
        else:
            # Handle different JSON structures
            continue
            
        # Construct full image path
        if local_path:
            img_path = os.path.join(root_path, local_path, filename)
        else:
            img_path = os.path.join(root_path, filename)
            
        if not os.path.exists(img_path):
            print(f"  Warning: Image not found: {img_path}")
            continue
            
        # Read EXIF
        try:
            with open(img_path, 'rb') as f:
                tags = exifread.process_file(f, details=False, stop_tag='DateTimeOriginal')
                timestamp = parse_exif_datetime(tags)
                
                if timestamp is not None:
                    timestamps[view_id] = timestamp
                    print(f"  View {view_id}: {datetime.fromtimestamp(timestamp).isoformat()} ({filename})")
                else:
                    print(f"  Warning: No timestamp found in EXIF for view {view_id}: {filename}")
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
        '--sfm_data',
        required=True,
        help='Input sfm_data.json file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output sfm_data.json file with constraints'
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
    
    # Load sfm_data
    print(f"Loading {args.sfm_data}...")
    with open(args.sfm_data, 'r') as f:
        sfm_data = json.load(f)
    
    # Get root path
    root_path = sfm_data.get('root_path', os.path.dirname(args.sfm_data))
    print(f"Root path: {root_path}")
    
    # Extract timestamps
    timestamps = extract_timestamps_from_images(sfm_data, root_path)
    
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
    
    # Save output
    print(f"\nSaving to {args.output}...")
    with open(args.output, 'w') as f:
        json.dump(sfm_data, f, indent=2)
    
    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
