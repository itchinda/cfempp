#!/bin/bash

# --- COLOR CONFIGURATION ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# --- 1. SECURITY CHECK ---
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${RED}[ERROR] Please run this script using 'source':${NC}"
    echo -e "      source scripts/setup_env.sh [release|debug]"
    exit 1
fi

echo -e "${BLUE}==> [CFEM++] Initializing Environment Setup...${NC}"

# --- 2. COMPILER & SYSTEM CHECKS ---
echo -e "==> Checking Compiler Capabilities..."
GCC_VERSION=$(g++ -dumpfullversion | cut -d'.' -f1)
if [ "$GCC_VERSION" -lt 10 ]; then
    echo -e "${RED}[ERROR] Your g++ version ($GCC_VERSION) is too old.${NC}"
    echo -e "CFEM++ requires at least g++ 10 for C++20 support."
    return 1
fi

if ! command -v pkg-config &> /dev/null; then
    echo -e "${BLUE}==> pkg-config not found. Installing...${NC}"
    sudo apt-get update && sudo apt-get install -y pkg-config
else
    echo -e "${GREEN}[OK] pkg-config is already installed.${NC}"
fi

echo -e "${BLUE}==> Installing system dependencies...${NC}"
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake python3 m4 wget git \
                        libgmp-dev libmpfr-dev libmpc-dev libzstd-dev \
                        libedit-dev libcurl4-openssl-dev llvm-dev libomp-dev zlib1g-dev

# --- 3. ARGUMENT PARSING & BUILD MODE ---
# Default to release if no argument provided
MODE_ARG=$(echo "${1:-release}" | tr '[:upper:]' '[:lower:]')

if [[ "$MODE_ARG" == "debug" ]]; then
    export CMAKE_BUILD_TYPE="Debug"
    export PETSC_ARCH="arch-linux-c-debug"
else
    export CMAKE_BUILD_TYPE="Release"
    export PETSC_ARCH="arch-linux-c-opt"
fi

echo -e "${BLUE}==> [CFEM++] Setting up environment for ${GREEN}${CMAKE_BUILD_TYPE}${BLUE} mode...${NC}"

# --- 4. ENVIRONMENT PATHS ---
USER_HOME=$(eval echo "~$USER")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"

# --- 5. SMART PETSC DETECTION ---
echo -e "==> Searching for PETSc installation ($PETSC_ARCH)..."

# Check if user already has a valid PETSC_DIR in their environment
if [ -n "$PETSC_DIR" ] && [ -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    echo -e "${GREEN}[OK] Using existing PETSC_DIR from environment: $PETSC_DIR${NC}"
else
    # Search the home directory up to 4 levels deep for the specific architecture configuration file
    SEARCH_PETSC=$(find "$USER_HOME" -maxdepth 4 -name "petscvariables" 2>/dev/null | grep "$PETSC_ARCH" | head -n 1)

    if [ -n "$SEARCH_PETSC" ]; then
        # Navigate up from .../lib/petsc/conf/petscvariables to the root PETSc directory
        export PETSC_DIR=$(cd "$(dirname "$SEARCH_PETSC")/../../../.." && pwd)
        echo -e "${GREEN}[OK] Automatically detected PETSC_DIR at: $PETSC_DIR${NC}"
    else
        export PETSC_DIR="$USER_HOME/projects/petsc"
        echo -e "${RED}[WARN] PETSc ($PETSC_ARCH) not found. Defaulting to: $PETSC_DIR${NC}"
    fi
fi

# --- 6. SMART SYMENGINE DETECTION ---
echo -e "==> Searching for SymEngine..."
IS_LOCAL_BUILD=false

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
    echo -e "${RED}[WARN] SymEngine not found. Defaulting to: $SYMENGINE_INSTALL_DIR${NC}"
fi

echo -e "Build Mode     : ${BLUE}$CMAKE_BUILD_TYPE${NC}"
echo -e "PETSC_DIR      : $PETSC_DIR"
echo -e "PETSC_ARCH     : $PETSC_ARCH"

# --- 7. AUTOMATED INSTALLATIONS (If missing) ---

# PETSc Installation
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    echo -e "${BLUE}==> PETSc architecture $PETSC_ARCH not detected. Cloning and building...${NC}"
    if [ ! -d "$PETSC_DIR" ]; then
        git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
    fi
    cd "$PETSC_DIR"
    
    # Enable debugging only if in Debug mode
    if [[ "$CMAKE_BUILD_TYPE" == "Debug" ]]; then
        DEBUG_FLAG="yes"
    else
        DEBUG_FLAG="no"
    fi

    ./configure --with-debugging=$DEBUG_FLAG \
                --download-fblaslapack \
                --download-mpich \
                --with-shared-libraries=1 \
                --PETSC_ARCH="$PETSC_ARCH"
                
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all
    cd -
fi

# SymEngine Installation
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

# --- 8. SESSION UPDATE ---
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# --- 9. CFEM++ CONFIGURATION ---
echo -e "${BLUE}==> Configuring CFEM++ Project with CMake...${NC}"
cd "$PROJECT_ROOT"
mkdir -p build && cd build
rm -f CMakeCache.txt # Clean cache to ensure clean mode switch

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SYMENGINE_CMAKE_PATH"

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY (${CMAKE_BUILD_TYPE})!"
echo -e "  Run 'cd build && make -j' to compile."
echo -e "=============================================${NC}"