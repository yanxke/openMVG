#!/usr/bin/env python3
"""
Analyze the heading angle distribution of generated pairs.
Shows histogram of heading differences between image pairs.
"""

import argparse
import json
import os
import math
from collections import defaultdict


def read_match_pairs_json(filepath):
    """Read match pairs from JSON file"""
    pairs = []
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        exit(1)

    with open(filepath, 'r') as f:
        data = json.load(f)

    if 'pairs' in data:
        for pair_data in data['pairs']:
            i = pair_data['i']
            j = pair_data['j']
            pairs.append((i, j))

    return pairs


def read_sfm_data(filepath):
    """Read sfm_data.json and return dict of {view_id: {heading: float}}"""
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        exit(1)

    with open(filepath, 'r') as f:
        data = json.load(f)

    view_data = {}

    # Parse views and their priors
    if 'views' in data:
        for view in data['views']:
            view_id = view['value']['ptr_wrapper']['data']['id_view']
            view_data[view_id] = {}

            view_info = view['value']['ptr_wrapper']['data']

            # Try multiple possible heading field names
            heading = None
            if 'gps_heading' in view_info and view_info.get('has_gps_heading', False):
                heading = view_info['gps_heading']
            elif 'heading' in view_info:
                heading = view_info['heading']
            elif 'center_pose_prior' in view_info:
                priors = view_info['center_pose_prior']
                if 'heading' in priors:
                    heading = priors['heading']

            if heading is not None:
                view_data[view_id]['heading'] = heading

    return view_data


def normalize_angle_diff(angle1, angle2):
    """
    Calculate the smallest angle difference between two headings (0-360).
    Returns value in range [0, 180].
    """
    diff = abs(angle1 - angle2)
    if diff > 180:
        diff = 360 - diff
    return diff


def create_histogram(angle_diffs, bucket_size=10):
    """Create histogram of angle differences"""
    histogram = defaultdict(int)

    for angle in angle_diffs:
        bucket = int(angle // bucket_size) * bucket_size
        histogram[bucket] += 1

    return histogram


def print_histogram(histogram, angle_diffs, bucket_size=10):
    """Print histogram as ASCII art"""
    if not histogram:
        print("No data to display")
        return

    total_matches = len(angle_diffs)
    max_count = max(histogram.values())
    bar_width = 60  # characters

    print("\n" + "="*80)
    print("HEADING ANGLE DIFFERENCE HISTOGRAM")
    print("="*80)
    print(f"\nTotal image pairs with heading data: {total_matches}")
    print(f"Bucket size: {bucket_size}°\n")

    # Print histogram
    for bucket in range(0, 180, bucket_size):
        count = histogram.get(bucket, 0)
        percentage = 100 * count / total_matches if total_matches > 0 else 0

        # Create bar
        bar_length = int((count / max_count) * bar_width) if max_count > 0 else 0
        bar = '█' * bar_length

        # Format label
        label = f"{bucket:3d}°-{bucket+bucket_size-1:3d}°"

        print(f"{label}: {bar} {count:5d} ({percentage:5.1f}%)")

    # Summary statistics
    print("\n" + "-"*80)
    print("SUMMARY STATISTICS")
    print("-"*80)

    if angle_diffs:
        angle_diffs_sorted = sorted(angle_diffs)
        mean_angle = sum(angle_diffs) / len(angle_diffs)
        median_angle = angle_diffs_sorted[len(angle_diffs_sorted) // 2]
        min_angle = min(angle_diffs)
        max_angle = max(angle_diffs)

        print(f"Mean:   {mean_angle:6.1f}°")
        print(f"Median: {median_angle:6.1f}°")
        print(f"Min:    {min_angle:6.1f}°")
        print(f"Max:    {max_angle:6.1f}°")

        # Cumulative percentages at key thresholds
        print("\n" + "-"*80)
        print("CUMULATIVE DISTRIBUTION")
        print("-"*80)

        thresholds = [30, 60, 90, 120, 150, 180]
        for threshold in thresholds:
            count_under = sum(1 for a in angle_diffs if a <= threshold)
            percentage = 100 * count_under / len(angle_diffs)
            print(f"Within {threshold:3d}°: {count_under:5d} ({percentage:5.1f}%)")

    print("\n" + "="*80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze heading angle distribution of verified matches.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s /data/task-2026-02-01/
  %(prog)s -p /data/task-2026-02-01/ -b 15
        """)

    parser.add_argument('project_folder',
                        help='Project folder containing openmvg/matches/')
    parser.add_argument('-b', '--bucket-size', type=int, default=10,
                        help='Histogram bucket size in degrees (default: 10)')

    args = parser.parse_args()

    # Construct file paths
    matches_json_file = os.path.join(args.project_folder, 'openmvg', 'matches', 'matches.e.json')
    sfm_data_file = os.path.join(args.project_folder, 'openmvg', 'matches', 'sfm_data.json')

    # Validate files exist
    print("Validating input files...")
    all_valid = True
    for filepath, name in [(matches_json_file, "matches.e.json"),
                            (sfm_data_file, "sfm_data.json")]:
        if not os.path.exists(filepath):
            print(f"Error: {name} not found: {filepath}")
            all_valid = False

    if not all_valid:
        print("\nPlease check the file paths and try again.")
        exit(1)

    # Read data
    print("Reading files...")
    pairs = read_match_pairs_json(matches_json_file)
    view_data = read_sfm_data(sfm_data_file)

    print(f"  Loaded {len(pairs):,} verified match image pairs")
    print(f"  Loaded {len(view_data):,} views from sfm_data.json")

    # Calculate angle differences
    angle_diffs = []
    pairs_without_heading = 0

    for i, j in pairs:
        # Check if both views have heading data
        if i in view_data and 'heading' in view_data[i] and \
           j in view_data and 'heading' in view_data[j]:

            heading_i = view_data[i]['heading']
            heading_j = view_data[j]['heading']

            angle_diff = normalize_angle_diff(heading_i, heading_j)
            angle_diffs.append(angle_diff)
        else:
            pairs_without_heading += 1

    if pairs_without_heading > 0:
        print(f"\nWarning: {pairs_without_heading} image pairs have missing heading data")

    # Create and print histogram
    if not angle_diffs:
        print("\nError: No image pairs with heading data found!")
        exit(1)

    histogram = create_histogram(angle_diffs, args.bucket_size)
    print_histogram(histogram, angle_diffs, args.bucket_size)


if __name__ == "__main__":
    main()
