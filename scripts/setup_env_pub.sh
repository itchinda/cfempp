#!/bin/bash
# ************************************************************************
# --- CFEM++ Library
# ---
# --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
# --- ALL RIGHTS RESERVED
# ---
# --- This software is protected by international copyright laws.
# --- Use, distribution, or modification of this software in any form,
# --- source or binary, for personal, academic, or non-commercial 
# --- purposes is permitted free of charge, provided that this 
# --- copyright notice and this permission notice appear in all 
# --- copies and supporting documentation.
# ---
# --- The authors and contributors provide this code "as is" without 
# --- any express or implied warranty. In no event shall they be 
# --- held liable for any damages arising from the use of this software.
# ************************************************************************

# ==============================================================================
# CFEM++ Environment Setup Script (Public Version)
# Verifies dependencies and configures the CMake build environment safely.
# Supports --bootstrap to automatically build missing dependencies.
# ==============================================================================

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Vérification que le script est bien sourcé
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${RED}[ERROR] Please run using 'source': source scripts/setup_env_pub.sh [release|debug] [--bootstrap]${NC}"
    exit 1
fi

echo -e "${BLUE}==> [CFEM++] Verifying Environment...${NC}"

# --- PARSE ARGUMENTS ---
MODE_ARG="release"
DO_BOOTSTRAP=false

for arg in "$@"; do
    case $(echo "$arg" | tr '[:upper:]' '[:lower:]') in
        debug|release)
            MODE_ARG="$arg"
            ;;
        --bootstrap)
            DO_BOOTSTRAP=true
            ;;
        *)
            echo -e "${YELLOW}[WARN] Unknown argument: $arg. Ignoring.${NC}"
            ;;
    esac
done

CMAKE_BUILD_TYPE=$([[ "$MODE_ARG" == "debug" ]] && echo "Debug" || echo "Release")
OS_PREFIX=$([[ "$(uname -s)" == "Darwin" ]] && echo "arch-darwin" || echo "arch-linux")
PETSC_ARCH=$([[ "$MODE_ARG" == "debug" ]] && echo "${OS_PREFIX}-c-debug" || echo "${OS_PREFIX}-c-opt")

# PATHS & HARDWARE
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"
SEARCH_BASE_DIRS=("$PARENT_DIR" "$HOME/software" "/opt" "/opt/homebrew" "/usr/local" "/usr/lib")
NCORES=$(command -v nproc &> /dev/null && nproc || sysctl -n hw.ncpu || echo 4)

# SYSTEM DEPENDENCIES INFO
if [ "$DO_BOOTSTRAP" = false ]; then
    echo -e "${YELLOW}[INFO] System Dependencies${NC}"
    echo -e "CFEM++ requires standard math and compilation libraries (GMP, MPFR, LLVM, etc.)."
    echo -e "Debian/Ubuntu: sudo apt-get install build-essential gfortran cmake git pkg-config libgmp-dev libmpfr-dev libmpc-dev libzstd-dev libcurl4-openssl-dev llvm-dev libomp-dev"
    echo -e "macOS (Homebrew): brew install gcc cmake git pkg-config gmp mpfr libmpc zstd curl llvm libomp"
    echo -e "------------------------------------------------------------"
fi

# SMART PETSC DETECTION & BOOTSTRAP
echo -e "${BLUE}==> Checking for PETSc ($PETSC_ARCH)...${NC}"
PETSC_FOUND=false

if [ -n "$PETSC_DIR" ] && { [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" ]] || [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib" ]] || [[ -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a" ]]; }; then
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
    if [ "$DO_BOOTSTRAP" = true ]; then
        echo -e "${YELLOW}[BOOTSTRAP] Building PETSc automatically...${NC}"
        INSTALL_DIR="$HOME/software"
        mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"
        
        [ ! -d "petsc" ] && git clone -b release https://gitlab.com/petsc/petsc.git petsc
        cd petsc
        
        # INTELLIGENT MPI DETECTION
        if command -v mpicc &> /dev/null && command -v mpicxx &> /dev/null; then
            echo -e "${GREEN}    -> System MPI detected, skipping MPICH download.${NC}"
            MPI_OPTION="--with-mpi=1"
            
            if command -v mpif90 &> /dev/null; then MPI_FC="mpif90"
            elif command -v mpifort &> /dev/null; then MPI_FC="mpifort"
            else MPI_FC="${FC:-gfortran}"; fi
            
            COMPILER_OPTS="--with-cc=mpicc --with-cxx=mpicxx --with-fc=$MPI_FC"
        else
            echo -e "${YELLOW}    -> No System MPI detected, forcing MPICH download.${NC}"
            MPI_OPTION="--download-mpich"
            COMPILER_OPTS="" # Let PETSc use system default compilers
        fi
        
        ./configure $COMPILER_OPTS \
                    --PETSC_ARCH="$PETSC_ARCH" \
                    --with-debugging=$([[ "$MODE_ARG" == "debug" ]] && echo "yes" || echo "no") \
                    --download-fblaslapack \
                    $MPI_OPTION \
                    --with-shared-libraries=1 || { echo -e "${RED}[ERROR] PETSc config failed${NC}"; cd "$PROJECT_ROOT"; return 1; }
        
        make PETSC_DIR=$(pwd) PETSC_ARCH="$PETSC_ARCH" -j$NCORES || { echo -e "${RED}[ERROR] PETSc build failed${NC}"; cd "$PROJECT_ROOT"; return 1; }
        
        export PETSC_DIR=$(pwd)
        cd "$PROJECT_ROOT"
    else
        echo -e "${RED}[ERROR] PETSc ($PETSC_ARCH) could not be found!${NC}"
        echo -e "${YELLOW}Run this script with '--bootstrap' to install it automatically, or do it manually:${NC}"
        echo -e "  git clone -b release https://gitlab.com/petsc/petsc.git"
        echo -e "  cd petsc && ./configure --PETSC_ARCH=$PETSC_ARCH \\
                                           --download-fblaslapack \\
                                           --download-mpich \\
                                           --with-shared-libraries=1"
        echo -e "  make PETSC_DIR=\$(pwd) PETSC_ARCH=$PETSC_ARCH all"
        return 1
    fi
fi
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# SMART SYMENGINE DETECTION & BOOTSTRAP
echo -e "${BLUE}==> Checking for SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
SYM_FOUND=false

ALL_FOUND_SYMS=$(find "${SEARCH_BASE_DIRS[@]}" -maxdepth 5 -name "SymEngineConfig.cmake" 2>/dev/null)

for FOUND_SYM in $ALL_FOUND_SYMS; do
    SYM_CONFIG_DIR=$(dirname "$FOUND_SYM")
    CACHE_FILE="$SYM_CONFIG_DIR/CMakeCache.txt"
    
    if [ -f "$CACHE_FILE" ]; then
        DETECTED_MODE=$(grep -i "^CMAKE_BUILD_TYPE" "$CACHE_FILE" | cut -d "=" -f 2 | tr '[:upper:]' | tr '[:lower:]' | tr -d ' ')
        if [ "$DETECTED_MODE" == "$MODE_ARG" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] SymEngine ($DETECTED_MODE) detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
    else
        if [ -f "$SYM_CONFIG_DIR/SymEngineTargets-${MODE_ARG}.cmake" ] || { [ "$MODE_ARG" == "release" ] && [ -f "$SYM_CONFIG_DIR/SymEngineTargets-noconfig.cmake" ]; }; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] System SymEngine detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
    fi
done

if [ "$SYM_FOUND" = false ]; then
    if [ "$DO_BOOTSTRAP" = true ]; then
        echo -e "${YELLOW}[BOOTSTRAP] Building SymEngine automatically...${NC}"
        INSTALL_DIR="$HOME/software"
        mkdir -p "$INSTALL_DIR" && cd "$INSTALL_DIR"
        
        [ ! -d "symengine" ] && git clone https://github.com/symengine/symengine symengine
        mkdir -p "symengine/build_${MODE_ARG}" && cd "symengine/build_${MODE_ARG}"
        
        cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
                 -DBUILD_SHARED_LIBS=ON \
                 -DWITH_SYMENGINE_RCP=OFF \
                 -DWITH_LLVM=ON \
                 -Dln_depend_mode=OFF || { echo -e "${RED}[ERROR] SymEngine config failed${NC}"; cd "$PROJECT_ROOT"; return 1; }
                 
        make -j$NCORES || { echo -e "${RED}[ERROR] SymEngine build failed${NC}"; cd "$PROJECT_ROOT"; return 1; }
        
        export SymEngine_DIR=$(pwd)
        cd "$PROJECT_ROOT"
    else
        echo -e "${RED}[ERROR] SymEngine could not be found for mode: $CMAKE_BUILD_TYPE!${NC}"
        echo -e "${YELLOW}Run this script with '--bootstrap' to install it automatically, or do it manually:${NC}"
        echo -e "  git clone https://github.com/symengine/symengine"
        echo -e "  mkdir symengine/build_${MODE_ARG} && cd symengine/build_${MODE_ARG}"
        echo -e "  cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DBUILD_SHARED_LIBS=ON -DWITH_LLVM=ON"
        echo -e "  make -j"
        return 1
    fi
fi

# --- CFEM++ CONFIGURATION ---
BUILD_DIR="build_${MODE_ARG}"
echo -e "${BLUE}==> Configuring CFEM++ Project with CMake in '$BUILD_DIR'...${NC}"
cd "$PROJECT_ROOT"
mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
rm -f CMakeCache.txt

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SymEngine_DIR" || { echo -e "${RED}[ERROR] CMake failed${NC}"; cd "$PROJECT_ROOT"; return 1; }

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY (${CMAKE_BUILD_TYPE})!"
echo -e "  Location: $PROJECT_ROOT/$BUILD_DIR"
echo -e "  Run 'cd $BUILD_DIR && make -j' to compile."
echo -e "=============================================${NC}"
cd "$PROJECT_ROOT"