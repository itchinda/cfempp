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

# Logic: Group all auto-installs in a single "software" directory
LOCAL_SOFTWARE_DIR="$USER_HOME/software"
mkdir -p "$LOCAL_SOFTWARE_DIR"

# --- 5. SMART PETSC DETECTION ---
echo -e "==> Searching for PETSc installation ($PETSC_ARCH)..."

if [ -n "$PETSC_DIR" ] && [ -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    echo -e "${GREEN}[OK] Using existing PETSC_DIR from environment: $PETSC_DIR${NC}"
else
    # Search Home for the variables file matching our current ARCH
    SEARCH_PETSC=$(find "$USER_HOME" -maxdepth 4 -name "petscvariables" 2>/dev/null | grep "$PETSC_ARCH" | head -n 1)

    if [ -n "$SEARCH_PETSC" ]; then
        # Backtrack to the PETSc root
        export PETSC_DIR=$(cd "$(dirname "$SEARCH_PETSC")/../../../.." && pwd)
        echo -e "${GREEN}[OK] Automatically detected PETSC_DIR at: $PETSC_DIR${NC}"
    else
        export PETSC_DIR="$LOCAL_SOFTWARE_DIR/petsc"
        echo -e "${RED}[WARN] PETSc ($PETSC_ARCH) not found. Fallback path: $PETSC_DIR${NC}"
    fi
fi

# --- 6. SMART SYMENGINE DETECTION ---
echo -e "==> Searching for SymEngine..."
IS_LOCAL_BUILD=false

# Look for the config file in sibling/parent directories
SEARCH_RESULTS=$(find "$PARENT_DIR" -maxdepth 4 -name "SymEngineConfig.cmake" 2>/dev/null)

if [ -n "$SEARCH_RESULTS" ]; then
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
    export SYMENGINE_INSTALL_DIR="$LOCAL_SOFTWARE_DIR/symengine_install"
    SYMENGINE_CMAKE_PATH="$SYMENGINE_INSTALL_DIR/lib/cmake/symengine"
    IS_LOCAL_BUILD=false
    echo -e "${RED}[WARN] SymEngine not found. Fallback path: $SYMENGINE_INSTALL_DIR${NC}"
fi

# Summary of chosen paths
echo -e "-------------------------------------------------------"
echo -e "Build Mode     : ${BLUE}$CMAKE_BUILD_TYPE${NC}"
echo -e "PETSC_DIR      : $PETSC_DIR"
echo -e "PETSC_ARCH     : $PETSC_ARCH"
echo -e "SymEngine Dir  : $SYMENGINE_INSTALL_DIR"
echo -e "-------------------------------------------------------"

# --- 7. AUTOMATED INSTALLATIONS ---

# PETSc Build Logic
if [ ! -d "$PETSC_DIR/$PETSC_ARCH" ]; then
    echo -e "${BLUE}==> PETSc architecture $PETSC_ARCH missing. Starting build...${NC}"
    if [ ! -d "$PETSC_DIR" ]; then
        git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
    fi
    cd "$PETSC_DIR"
    
    DEBUG_FLAG=$(if [[ "$CMAKE_BUILD_TYPE" == "Debug" ]]; then echo "yes"; else echo "no"; fi)

    ./configure --with-debugging=$DEBUG_FLAG \
                --download-fblaslapack \
                --download-mpich \
                --with-shared-libraries=1 \
                --PETSC_ARCH="$PETSC_ARCH"
                
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all
    cd -
fi

# SymEngine Build Logic
if [ ! -d "$SYMENGINE_INSTALL_DIR" ]; then
    echo -e "${BLUE}==> SymEngine not found. Installing into $SYMENGINE_INSTALL_DIR...${NC}"
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
rm -f CMakeCache.txt # Important for switching Release/Debug

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SYMENGINE_CMAKE_PATH"

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY (${CMAKE_BUILD_TYPE})!"
echo -e "  Location: $PROJECT_ROOT/build"
echo -e "  Run 'make -j' to compile."
echo -e "=============================================${NC}"