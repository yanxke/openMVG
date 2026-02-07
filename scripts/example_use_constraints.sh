#!/bin/bash
# Quick start guide for using time-based translation constraints

# Step 1: Extract EXIF timestamps and generate constraints
# This creates constraints based on 1.0 m/s maximum speed
python /home/yke/src/openMVG/scripts/populate_translation_constraints.py \
  --sfm_data matches/sfm_data.json \
  --output matches/sfm_data_with_constraints.json \
  --max_velocity 1.0 \
  --time_window 10.0

# Step 2: Run the global SfM pipeline with the constrained data
# Use -t 2 for SOFTL1 translation averaging (required for constraints)
openMVG_main_SfM \
  -i matches/sfm_data_with_constraints.json \
  -m matches \
  -o reconstruction \
  -t 2

# The solver will automatically:
# - Detect the max_translation_distance constraints
# - Remap them to the solver's coordinate system
# - Add penalty terms for violations
# - Log how many constraints were applied

# Expected output:
# "Using SOFTL1 solver with XXX max distance constraints."
# "Successfully added XXX distance constraints to the problem."
