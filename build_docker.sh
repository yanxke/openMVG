#!/bin/bash
set -e

# Check if repository is clean (no uncommitted changes)
if [ -z "$(git status --porcelain 2>/dev/null)" ]; then
  # Repository is clean - capture git short hash
  GIT_HASH=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")
  echo "Building Docker image with git hash: $GIT_HASH (repo is clean)"
  BUILD_ARGS="--build-arg GIT_SHORT_HASH=$GIT_HASH"
else
  # Repository has uncommitted changes - don't include git hash
  echo "Building Docker image without git hash (repo has uncommitted changes)"
  BUILD_ARGS=""
fi

# Build Docker image
docker build \
  $BUILD_ARGS \
  -t openmvg \
  .

echo ""
echo "Docker image built successfully!"
echo "  Tagged as: openmvg"
echo ""
echo "Run with: docker run openmvg"
