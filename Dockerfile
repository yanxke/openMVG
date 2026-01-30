# Use Ubuntu 22.04 (will be supported until April 2027)
FROM ubuntu:jammy

# Add openMVG binaries to path
ENV PATH=$PATH:/opt/openMVG_Build/install/bin

# Get dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata; \
  apt-get install -y \
  cmake \
  build-essential \
  ccache \
  graphviz \
  git \
  coinor-libclp-dev \
  libceres-dev \
  libflann-dev \
  libjpeg-dev \
  liblemon-dev \
  libpng-dev \
  libtiff-dev \
  libqt5svg5-dev \
  qttools5-dev \
  git \
  python3-minimal \
  python3-pip; \
  apt-get autoclean && apt-get clean

# Python deps for pose prior helper
RUN python3 -m pip install --no-cache-dir piexif

# Clone the openvMVG repo
ADD . /opt/openMVG

# Build
RUN --mount=type=cache,target=/ccache \
  mkdir -p /opt/openMVG_Build; \
  cd /opt/openMVG_Build; \
  export CCACHE_DIR=/ccache; \
  export CCACHE_BASEDIR=/opt/openMVG; \
  cmake -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX="/opt/openMVG_Build/install" \
    -DOpenMVG_BUILD_TESTS=OFF \
    -DOpenMVG_BUILD_EXAMPLES=OFF \
    -DOpenMVG_BUILD_GUI_SOFTWARES=OFF \
    -DOpenMVG_BUILD_DOCS=OFF \
    -DCMAKE_C_COMPILER_LAUNCHER=ccache \
    -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DCOINUTILS_INCLUDE_DIR_HINTS=/usr/include \
    -DLEMON_INCLUDE_DIR_HINTS=/usr/include/lemon \
    -DCLP_INCLUDE_DIR_HINTS=/usr/include \
    -DOSI_INCLUDE_DIR_HINTS=/usr/include \
    ../openMVG/src; \
    make -j 8;

RUN --mount=type=cache,target=/ccache \
  cd /opt/openMVG_Build && make install;
