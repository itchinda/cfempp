#!/bin/bash
# ==============================================================================
# CFEM++ Environment Setup Script (Public Version)
# Verifies dependencies and configures the CMake build environment safely.
# ==============================================================================

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${RED}[ERROR] Please run using 'source': source scripts/setup_env_pub.sh [release|debug]${NC}"
    exit 1
fi

echo -e "${BLUE}==> [CFEM++] Verifying Environment...${NC}"

# --- 1. SYSTEM DEPENDENCIES INFO ---
echo -e "${YELLOW}[INFO] System Dependencies${NC}"
echo -e "CFEM++ requires standard math and compilation libraries (GMP, MPFR, LLVM, etc.)."
echo -e "For Debian/Ubuntu users, you can install them via:"
echo -e "  sudo apt-get install build-essential gfortran cmake git pkg-config \\"
echo -e "       libgmp-dev libmpfr-dev libmpc-dev libzstd-dev \\"
echo -e "       libcurl4-openssl-dev llvm-dev libomp-dev"
echo -e "------------------------------------------------------------"

# --- MODE ---
MODE_ARG=$(echo "${1:-release}" | tr '[:upper:]' '[:lower:]')
CMAKE_BUILD_TYPE=$([[ "$MODE_ARG" == "debug" ]] && echo "Debug" || echo "Release")
PETSC_ARCH=$([[ "$MODE_ARG" == "debug" ]] && echo "arch-linux-c-debug" || echo "arch-linux-c-opt")

# --- PATHS ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"
SEARCH_BASE_DIRS=("$PARENT_DIR" "$HOME/software" "/opt" "/usr/local" "/usr/lib")

# --- SMART PETSC DETECTION ---
echo -e "${BLUE}==> Checking for PETSc ($PETSC_ARCH)...${NC}"
PETSC_FOUND=false

if [ -n "$PETSC_DIR" ] && { [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" ]] || [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a" ]]; }; then
    PETSC_FOUND=true
else
    for base_dir in "${SEARCH_BASE_DIRS[@]}"; do
        if [ -d "$base_dir" ]; then
            FOUND_PATHS=$(find "$base_dir" -maxdepth 5 -name "petscvariables" 2>/dev/null | grep "$PETSC_ARCH" | head -n 1)
            if [ -n "$FOUND_PATHS" ]; then
                export PETSC_DIR=$(cd "$(dirname "$FOUND_PATHS")/../../../.." && pwd)
                PETSC_FOUND=true; break
            fi
        fi
    done
fi

if [ "$PETSC_FOUND" = false ]; then
    echo -e "${RED}[ERROR] PETSc ($PETSC_ARCH) could not be found!${NC}"
    echo -e "${YELLOW}Please install PETSc for your target mode (${CMAKE_BUILD_TYPE}):${NC}"
    echo -e "  git clone -b release https://gitlab.com/petsc/petsc.git"
    echo -e "  cd petsc"
    echo -e "  ./configure --PETSC_ARCH=$PETSC_ARCH \\"
    echo -e "              --with-debugging=$([[ "$MODE_ARG" == "debug" ]] && echo "yes" || echo "no") \\"
    echo -e "              --download-fblaslapack \\"
    echo -e "              --download-mpich \\"
    echo -e "              --with-shared-libraries=1"
    echo -e "  make PETSC_DIR=\$(pwd) PETSC_ARCH=$PETSC_ARCH all"
    echo -e "\nThen run: export PETSC_DIR=/path/to/petsc before running this script again."
    return 1
fi
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# --- SMART SYMENGINE DETECTION (CACHE BASED) ---
echo -e "${BLUE}==> Checking for SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
SYM_FOUND=false

ALL_FOUND_SYMS=$(find "${SEARCH_BASE_DIRS[@]}" -maxdepth 5 -name "SymEngineConfig.cmake" 2>/dev/null)

for FOUND_SYM in $ALL_FOUND_SYMS; do
    SYM_CONFIG_DIR=$(dirname "$FOUND_SYM")
    CACHE_FILE="$SYM_CONFIG_DIR/CMakeCache.txt"
    
    # Cas A : Build local (CMakeCache.txt existe)
    if [ -f "$CACHE_FILE" ]; then
        DETECTED_MODE=$(grep -i "^CMAKE_BUILD_TYPE" "$CACHE_FILE" | cut -d "=" -f 2 | tr '[:upper:]' '[:lower:]' | tr -d ' ')
        if [ "$DETECTED_MODE" == "$MODE_ARG" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] SymEngine ($DETECTED_MODE) detected in build tree at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
        
    # Cas B : Installation système globale (make install)
    else
        if [ -f "$SYM_CONFIG_DIR/SymEngineTargets-${MODE_ARG}.cmake" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] System SymEngine (matches $MODE_ARG) detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
            
        elif [ "$MODE_ARG" == "release" ] && [ -f "$SYM_CONFIG_DIR/SymEngineTargets-noconfig.cmake" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] System SymEngine (noconfig -> Release) detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
    fi
done

if [ "$SYM_FOUND" = false ]; then
    echo -e "${RED}[ERROR] SymEngine could not be found for mode: $CMAKE_BUILD_TYPE!${NC}"
    echo -e "${YELLOW}Please build SymEngine with LLVM support:${NC}"
    echo -e "  git clone https://github.com/symengine/symengine"
    echo -e "  mkdir symengine/build_${MODE_ARG} && cd symengine/build_${MODE_ARG}"
    echo -e "  cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \\"
    echo -e "           -DBUILD_SHARED_LIBS=ON \\"
    echo -e "           -DWITH_SYMENGINE_RCP=OFF \\"
    echo -e "           -DWITH_LLVM=ON \\"
    echo -e "           -Dln_depend_mode=OFF"
    echo -e "  make -j"
    echo -e "\nThen run: export SymEngine_DIR=/path/to/symengine/build_${MODE_ARG} before running this script again."
    return 1
fi

# --- CFEM++ CONFIGURATION ---
echo -e "${BLUE}==> Configuring CFEM++ Project with CMake...${NC}"
cd "$PROJECT_ROOT"
mkdir -p build && cd build
rm -f CMakeCache.txt

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SymEngine_DIR" || { echo -e "${RED}[ERROR] CMake failed${NC}"; return 1; }

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY (${CMAKE_BUILD_TYPE})!"
echo -e "  Location: $PROJECT_ROOT/build"
echo -e "  Run 'cd build && make -j' to compile."
echo -e "=============================================${NC}"
