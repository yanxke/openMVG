#!/bin/bash
set -e

# Capture git short hash
GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

echo "Building Docker image with git hash: $GIT_HASH"

# Build Docker image with git hash as build argument
docker build \
  --build-arg GIT_SHORT_HASH="$GIT_HASH" \
  -t openmvg:latest \
  .

echo ""
echo "Docker image built successfully!"
echo "  Tagged as: openmvg:latest"
echo ""
echo "Run with: docker run openmvg:latest"
