#!/bin/bash

# --- COLOR CONFIGURATION ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# --- SECURITY CHECK ---
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${RED}[ERROR] Please run this script using 'source':${NC}"
    echo -e "      source scripts/setup_env.sh [release|debug]"
    exit 1
fi

echo -e "${BLUE}==> [CFEM++] Initializing Environment Setup...${NC}"

# --- COMPILER & SYSTEM CHECKS ---
echo -e "==> Checking Compiler Capabilities..."
GCC_VERSION=$(g++ -dumpfullversion | cut -d'.' -f1)
if [ "$GCC_VERSION" -lt 10 ]; then
    echo -e "${RED}[ERROR] Your g++ version ($GCC_VERSION) is too old.${NC}"
    echo -e "CFEM++ requires at least g++ 10 for C++20 support."
    exit 1
fi

if ! command -v pkg-config &> /dev/null; then
    echo -e "${BLUE}==> pkg-config not found. Installing...${NC}"
    sudo apt-get update && sudo apt-get install -y pkg-config
else
    echo -e "${GREEN}[OK] pkg-config is already installed.${NC}"
fi

echo -e "${BLUE}==> Installing system dependencies...${NC}"
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake m4 wget git \
                        libgmp-dev libmpfr-dev libmpc-dev libzstd-dev \
                        libedit-dev libcurl4-openssl-dev llvm-dev libomp-dev zlib1g-dev

# --- ARGUMENT PARSING & BUILD MODE ---
# Default to release if no argument provided
MODE_ARG=$(echo "${1:-release}" | tr '[:upper:]' '[:lower:]')

if [[ "$MODE_ARG" == "debug" ]]; then
    export CMAKE_BUILD_TYPE="Debug"
    DEFAULT_PETSC_ARCH="arch-linux-c-debug"
    LOOK_FOR="build_debug"
else
    export CMAKE_BUILD_TYPE="Release"
    DEFAULT_PETSC_ARCH="arch-linux-c-opt"
    LOOK_FOR="build_release"
fi

# --- ENVIRONMENT PATHS ---
USER_HOME=$(eval echo "~$USER")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"

# --- SMART SYMENGINE DETECTION ---
IS_LOCAL_BUILD=false
echo -e "==> Searching for SymEngine..."

# Search for the config file fingerprint
SEARCH_RESULTS=$(find "$PARENT_DIR" -maxdepth 4 -name "SymEngineConfig.cmake" 2>/dev/null)

if [ -n "$SEARCH_RESULTS" ]; then
    # Prioritize match for current mode
    MATCHING_PATH=$(echo "$SEARCH_RESULTS" | grep -i "$MODE_ARG" | head -n 1)
    FINAL_CONFIG=$(if [ -n "$MATCHING_PATH" ]; then echo "$MATCHING_PATH"; else echo "$SEARCH_RESULTS" | head -n 1; fi)
    
    CONFIG_DIR=$(dirname "$FINAL_CONFIG")
    if [[ "$CONFIG_DIR" == *"lib/cmake"* ]]; then
        export SYMENGINE_INSTALL_DIR=$(cd "$CONFIG_DIR/../../.." && pwd)
        SYMENGINE_CMAKE_PATH="$CONFIG_DIR"
        IS_LOCAL_BUILD=false
    else
        export SYMENGINE_INSTALL_DIR=$(cd "$CONFIG_DIR" && pwd)
        SYMENGINE_CMAKE_PATH="$SYMENGINE_INSTALL_DIR"
        IS_LOCAL_BUILD=true
    fi
    echo -e "${GREEN}[OK] SymEngine found at: $SYMENGINE_INSTALL_DIR${NC}"
else
    export SYMENGINE_INSTALL_DIR="$USER_HOME/symengine_install"
    SYMENGINE_CMAKE_PATH="$SYMENGINE_INSTALL_DIR/lib/cmake/symengine"
    IS_LOCAL_BUILD=false
fi

# --- PETSC SETTINGS ---
export PETSC_ARCH=${PETSC_ARCH:-$DEFAULT_PETSC_ARCH}
export PETSC_DIR=${PETSC_DIR:-$USER_HOME/projects/petsc}

echo -e "Build Mode     : ${BLUE}$CMAKE_BUILD_TYPE${NC}"
echo -e "PETSC_DIR      : $PETSC_DIR"
echo -e "PETSC_ARCH     : $PETSC_ARCH"

# --- 3. PETSC INSTALLATION ---
if [ ! -d "$PETSC_DIR" ]; then
    echo -e "${BLUE}==> PETSc not detected. Cloning and building...${NC}"
    git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
    cd "$PETSC_DIR"
    ./configure --with-debugging=no --download-fblaslapack --download-mpich --with-shared-libraries=1 --PETSC_ARCH="$PETSC_ARCH"
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all
    cd -
else
    echo -e "${GREEN}[OK] PETSc found.${NC}"
fi

# --- 4. SYMENGINE INSTALLATION ---
if [ ! -d "$SYMENGINE_INSTALL_DIR" ]; then
    echo -e "${BLUE}==> SymEngine not detected. Cloning and building...${NC}"
    git clone https://github.com/symengine/symengine ~/symengine_source
    mkdir -p ~/symengine_source/build && cd ~/symengine_source/build
    cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
             -DCMAKE_INSTALL_PREFIX="$SYMENGINE_INSTALL_DIR" \
             -DWITH_SYMENGINE_RCP=ON -DWITH_LLVM=ON -Dln_depend_mode=OFF
    make -j$(nproc) install
    cd - && rm -rf ~/symengine_source
    IS_LOCAL_BUILD=false
fi

# --- 5. SESSION UPDATE ---
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# --- 6. CFEM++ CONFIGURATION ---
echo -e "${BLUE}==> Configuring CFEM++ Project with CMake...${NC}"
mkdir -p build && cd build
rm -f CMakeCache.txt # Clean cache to ensure mode switch

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SYMENGINE_CMAKE_PATH"

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY (${CMAKE_BUILD_TYPE})!"
echo -e "  Run 'cd build && make -j' to compile."
echo -e "=============================================${NC}"