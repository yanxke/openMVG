#!/usr/bin/env python3
"""
Compare two pair generation methods to evaluate filtering effectiveness.
Shows how many pairs were filtered out and whether they would have passed geometric verification.
"""

import argparse
import struct
import os

def read_pairs_bin(filepath):
    """Read pairs from a .bin file and return as set of tuples (i, j) where i < j"""
    pairs = set()
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        exit(1)

    with open(filepath, 'rb') as f:
        # Read number of pairs
        num_pairs_bytes = f.read(8)
        if len(num_pairs_bytes) < 8:
            return pairs
        num_pairs = struct.unpack('Q', num_pairs_bytes)[0]

        # Read each pair
        for _ in range(num_pairs):
            pair_bytes = f.read(8)
            if len(pair_bytes) < 8:
                break
            i, j = struct.unpack('II', pair_bytes)
            # Normalize so smaller index is first
            pairs.add((min(i, j), max(i, j)))

    return pairs

def read_matches_bin(filepath):
    """Read geometrically verified matches and return as set of tuples (i, j) where i < j"""
    matches = set()
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found")
        exit(1)

    with open(filepath, 'rb') as f:
        # Read number of pairs
        num_pairs_bytes = f.read(8)
        if len(num_pairs_bytes) < 8:
            return matches
        num_pairs = struct.unpack('Q', num_pairs_bytes)[0]

        # Read each match
        for _ in range(num_pairs):
            # Read pair indices
            pair_bytes = f.read(8)
            if len(pair_bytes) < 8:
                break
            i, j = struct.unpack('II', pair_bytes)

            # Read number of matches for this pair
            num_matches_bytes = f.read(8)
            if len(num_matches_bytes) < 8:
                break
            num_matches = struct.unpack('Q', num_matches_bytes)[0]

            # Skip the actual match data (each match is 2 uint32s)
            f.seek(num_matches * 8, 1)

            # Normalize so smaller index is first
            matches.add((min(i, j), max(i, j)))

    return matches

def print_statistics(old_pairs, new_pairs, verified_matches):
    """Print comparison statistics"""

    print("\n" + "="*80)
    print("PAIR GENERATION COMPARISON")
    print("="*80)

    # Basic counts
    print(f"\nTotal pairs:")
    print(f"  Old method:                {len(old_pairs):>8,}")
    print(f"  New method:                {len(new_pairs):>8,}")
    print(f"  Reduction:                 {len(old_pairs) - len(new_pairs):>8,} ({100*(len(old_pairs)-len(new_pairs))/len(old_pairs):.1f}%)")

    print(f"\nGeometrically verified:")
    print(f"  Total verified matches:    {len(verified_matches):>8,}")

    # Set operations
    only_in_old = old_pairs - new_pairs
    only_in_new = new_pairs - old_pairs
    in_both = old_pairs & new_pairs

    print(f"\nSet comparison:")
    print(f"  In both methods:           {len(in_both):>8,}")
    print(f"  Only in old:               {len(only_in_old):>8,}")
    print(f"  Only in new:               {len(only_in_new):>8,}")

    # Match verification rates
    old_verified = old_pairs & verified_matches
    new_verified = new_pairs & verified_matches

    print(f"\n" + "-"*80)
    print("VERIFICATION RATES")
    print("-"*80)

    old_rate = 100 * len(old_verified) / len(old_pairs) if old_pairs else 0
    new_rate = 100 * len(new_verified) / len(new_pairs) if new_pairs else 0

    print(f"\nOld method:")
    print(f"  Pairs that verified:       {len(old_verified):>8,} / {len(old_pairs):>8,} ({old_rate:.1f}%)")
    print(f"  Pairs that failed:         {len(old_pairs) - len(old_verified):>8,} / {len(old_pairs):>8,} ({100-old_rate:.1f}%)")

    print(f"\nNew method:")
    print(f"  Pairs that verified:       {len(new_verified):>8,} / {len(new_pairs):>8,} ({new_rate:.1f}%)")
    print(f"  Pairs that failed:         {len(new_pairs) - len(new_verified):>8,} / {len(new_pairs):>8,} ({100-new_rate:.1f}%)")

    # Critical analysis: what did we miss?
    filtered_out = only_in_old  # Pairs in old but not new
    filtered_out_verified = filtered_out & verified_matches  # Pairs we filtered out that would have verified
    filtered_out_failed = filtered_out - verified_matches    # Pairs we filtered out that would have failed

    print(f"\n" + "-"*80)
    print("FILTERING EFFECTIVENESS (pairs removed by new method)")
    print("-"*80)

    print(f"\nPairs filtered out by new method: {len(filtered_out):>8,}")

    if filtered_out:
        miss_rate = 100 * len(filtered_out_verified) / len(filtered_out)
        save_rate = 100 * len(filtered_out_failed) / len(filtered_out)

        print(f"  Would have VERIFIED (FALSE NEGATIVES - BAD):  {len(filtered_out_verified):>8,} ({miss_rate:.1f}%)")
        print(f"  Would have FAILED (TRUE NEGATIVES - GOOD):    {len(filtered_out_failed):>8,} ({save_rate:.1f}%)")

    # Recovery rate: how many verified matches did we keep?
    total_possible_matches = verified_matches & old_pairs  # Matches that old method could find
    recovered_matches = verified_matches & new_pairs        # Matches that new method found
    missed_matches = verified_matches - new_pairs           # Verified matches we missed

    print(f"\n" + "-"*80)
    print("MATCH RECOVERY (verified matches found)")
    print("-"*80)

    recovery_rate = 100 * len(recovered_matches) / len(total_possible_matches) if total_possible_matches else 0

    print(f"\nVerified matches found by old method:  {len(total_possible_matches):>8,}")
    print(f"Verified matches found by new method:  {len(recovered_matches):>8,} ({recovery_rate:.1f}%)")
    print(f"Verified matches MISSED by new method: {len(missed_matches):>8,} ({100-recovery_rate:.1f}%)")

    # Efficiency metrics
    print(f"\n" + "-"*80)
    print("EFFICIENCY SUMMARY")
    print("-"*80)

    old_efficiency = len(old_verified) / len(old_pairs) if old_pairs else 0
    new_efficiency = len(new_verified) / len(new_pairs) if new_pairs else 0

    print(f"\nPairs per verified match:")
    print(f"  Old method: {len(old_pairs) / len(old_verified):.2f} pairs/match" if old_verified else "  Old method: N/A")
    print(f"  New method: {len(new_pairs) / len(new_verified):.2f} pairs/match" if new_verified else "  New method: N/A")

    print(f"\nEfficiency improvement: {new_efficiency / old_efficiency:.2f}x" if old_efficiency > 0 else "N/A")
    print(f"Computational savings: {100 * (len(old_pairs) - len(new_pairs)) / len(old_pairs):.1f}% fewer pairs to process")

    # Verdict
    print(f"\n" + "="*80)
    print("VERDICT")
    print("="*80)

    if len(missed_matches) == 0:
        print("\n✓ PERFECT: New method found ALL verified matches with fewer pairs!")
    elif len(missed_matches) < 0.01 * len(total_possible_matches):
        print(f"\n✓ EXCELLENT: New method found {recovery_rate:.1f}% of verified matches")
        print(f"  Only {len(missed_matches)} matches missed (less than 1%)")
    elif len(missed_matches) < 0.05 * len(total_possible_matches):
        print(f"\n~ GOOD: New method found {recovery_rate:.1f}% of verified matches")
        print(f"  {len(missed_matches)} matches missed ({100-recovery_rate:.1f}%)")
    else:
        print(f"\n✗ NEEDS TUNING: New method missed {100-recovery_rate:.1f}% of verified matches")
        print(f"  Consider relaxing the heading threshold")

    if new_efficiency > old_efficiency * 1.2:
        print(f"✓ New method is more efficient ({new_rate:.1f}% vs {old_rate:.1f}% verification rate)")

    print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare two pair generation methods to evaluate filtering effectiveness.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare EXHAUSTIVE vs COMPASS pair generation
  %(prog)s \\
    -op /data/openmvg/matches/pairs_exhaustive.bin \\
    -np /data/openmvg/matches/pairs_compass.bin \\
    -m /data/openmvg/matches/matches.e.bin

  # Compare different heading thresholds
  %(prog)s \\
    --old-pairs /data/openmvg/matches/pairs_h90.bin \\
    --new-pairs /data/openmvg/matches/pairs_h120.bin \\
    --matches /data/openmvg/matches/matches.e.bin
        """)

    parser.add_argument('-op', '--old-pairs', required=True,
                        help='Old/baseline pairs.bin file (e.g., EXHAUSTIVE or higher threshold)')
    parser.add_argument('-np', '--new-pairs', required=True,
                        help='New/filtered pairs.bin file (e.g., COMPASS or lower threshold)')
    parser.add_argument('-m', '--matches', required=True,
                        help='Geometrically verified matches.e.bin file')

    args = parser.parse_args()

    # Validate all files exist before processing
    print("Validating input files...")
    all_valid = True
    for filepath, name in [(args.old_pairs, "old-pairs"),
                            (args.new_pairs, "new-pairs"),
                            (args.matches, "matches")]:
        if not os.path.exists(filepath):
            print(f"Error: {name} file not found: {filepath}")
            all_valid = False

    if not all_valid:
        print("\nPlease check the file paths and try again.")
        exit(1)

    print("Reading files...")
    old_pairs = read_pairs_bin(args.old_pairs)
    new_pairs = read_pairs_bin(args.new_pairs)
    verified_matches = read_matches_bin(args.matches)

    print(f"  Old pairs: {len(old_pairs):,} pairs")
    print(f"  New pairs: {len(new_pairs):,} pairs")
    print(f"  Verified: {len(verified_matches):,} matches")

    print_statistics(old_pairs, new_pairs, verified_matches)
