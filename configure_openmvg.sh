#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${OPENMVG_BUILD_DIR:-${ROOT_DIR}/openMVG_Build}"
INSTALL_DIR="${OPENMVG_INSTALL_PREFIX:-${BUILD_DIR}/install}"
BUILD_TYPE="${CMAKE_BUILD_TYPE:-RELEASE}"

mkdir -p "${BUILD_DIR}"

cmake_args=(
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}"
  -DOpenMVG_BUILD_TESTS=OFF
  -DOpenMVG_BUILD_EXAMPLES=OFF
  -DOpenMVG_BUILD_GUI_SOFTWARES=OFF
  -DOpenMVG_BUILD_DOCS=OFF
)

if command -v ccache >/dev/null 2>&1; then
  cmake_args+=(
    -DCMAKE_C_COMPILER_LAUNCHER=ccache
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache
  )
else
  echo "ccache not found; configuring without compiler launcher."
fi

echo "Configuring openMVG"
echo "  Source : ${ROOT_DIR}/src"
echo "  Build  : ${BUILD_DIR}"
echo "  Install: ${INSTALL_DIR}"

cmake -S "${ROOT_DIR}/src" -B "${BUILD_DIR}" "${cmake_args[@]}" "$@"

