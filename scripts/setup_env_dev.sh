#!/bin/bash
#************************************************************************
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
#************************************************************************

# ==============================================================================
# CFEM++ Environment Setup Script 
# Verifies dependencies and configures the CMake build environment safely.
# ==============================================================================

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}==> CFEM++ DEPENDENCIES BOOTSTRAP${NC}"

# ------------------------------------------------------------
# MODE & ARGUMENTS
# ------------------------------------------------------------
MODE_ARG=$(echo "${1:-release}" | tr '[:upper:]' '[:lower:]')
FORCE_UPDATE=false
if [[ "$2" == "--force" ]]; then
    FORCE_UPDATE=true
fi

CMAKE_BUILD_TYPE=$([[ "$MODE_ARG" == "debug" ]] && echo "Debug" || echo "Release")

# OS detection for PETSc the target name
OS_PREFIX=$([[ "$(uname -s)" == "Darwin" ]] && echo "arch-darwin" || echo "arch-linux")
PETSC_ARCH=$([[ "$MODE_ARG" == "debug" ]] && echo "${OS_PREFIX}-c-debug" || echo "${OS_PREFIX}-c-opt")

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"
INSTALL_DIR="$HOME/software"
mkdir -p "$INSTALL_DIR"

# ------------------------------------------------------------
# CORE DETECTION
# ------------------------------------------------------------
NCORES=$(command -v nproc &> /dev/null && nproc || sysctl -n hw.ncpu || echo 4)

# ------------------------------------------------------------
# PLATFORM DETECTION
# ------------------------------------------------------------
detect_package_manager() {
    if command -v apt-get &> /dev/null; then echo "apt"
    elif command -v dnf &> /dev/null; then echo "dnf"
    elif command -v pacman &> /dev/null; then echo "pacman"
    elif command -v brew &> /dev/null; then echo "brew"
    else echo "unknown"
    fi
}

# ------------------------------------------------------------
# SYSTEM DEPS & CHECKER
# ------------------------------------------------------------
PKG_MANAGER=$(detect_package_manager)
echo -e "${BLUE}==> Detected package manager: $PKG_MANAGER${NC}"
echo -e "${BLUE}==> Verifying system dependencies${NC}"

check_pkg() {
    case "$PKG_MANAGER" in
        apt) dpkg -s "$1" >/dev/null 2>&1 ;;
        dnf) rpm -q "$1" >/dev/null 2>&1 ;;
        pacman) pacman -Qi "$1" >/dev/null 2>&1 ;;
        brew) brew list "$1" >/dev/null 2>&1 ;;
        *) return 1 ;;
    esac
}

case "$PKG_MANAGER" in
    apt) PKGS=(build-essential gfortran cmake git pkg-config libgmp-dev libmpfr-dev libmpc-dev libzstd-dev libcurl4-openssl-dev llvm-dev libomp-dev) ;;
    dnf) PKGS=(gcc gcc-c++ gcc-gfortran cmake git pkgconf-pkg-config gmp-devel mpfr-devel libmpc-devel zstd-devel libcurl-devel llvm-devel libomp-devel) ;;
    pacman) PKGS=(base-devel gcc-fortran cmake git pkgconf gmp mpfr libmpc zstd curl llvm openmp) ;;
    brew) PKGS=(gcc cmake git pkg-config gmp mpfr libmpc zstd curl llvm libomp) ;;
    *) PKGS=() ;;
esac

MISSING_PKGS=()
for pkg in "${PKGS[@]}"; do
    if ! check_pkg "$pkg"; then
        MISSING_PKGS+=("$pkg")
    fi
done

# For MacOs, check 'gcc', not 'g++' for the compiler tool
COMPILER_TOOL=$([[ "$PKG_MANAGER" == "brew" ]] && echo "gcc" || echo "g++")
REQUIRED_TOOLS=($COMPILER_TOOL cmake git pkg-config make)

MISSING_TOOLS=0
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v $tool &> /dev/null; then
        MISSING_TOOLS=1
    fi
done

if [ ${#MISSING_PKGS[@]} -gt 0 ] || [ "$MISSING_TOOLS" -eq 1 ]; then
    echo -e "${BLUE}==> Installing missing dependencies...${NC}"
    case "$PKG_MANAGER" in
        apt) sudo apt-get update && sudo apt-get install -y "${MISSING_PKGS[@]}" "${REQUIRED_TOOLS[@]}" || exit 1 ;;
        dnf) sudo dnf install -y "${MISSING_PKGS[@]}" "${REQUIRED_TOOLS[@]}" || exit 1 ;;
        pacman) sudo pacman -Sy --noconfirm "${MISSING_PKGS[@]}" "${REQUIRED_TOOLS[@]}" || exit 1 ;;
        brew) brew install "${MISSING_PKGS[@]}" || exit 1 ;;
        *) echo -e "${YELLOW}[WARN] Please install dependencies manually.${NC}" ;;
    esac
else
    echo -e "${GREEN}==> All required tools and libraries found.${NC}"
fi

# ------------------------------------------------------------
# COMPILER SELECTION & SAFETY CHECKS
# ------------------------------------------------------------
if [[ "$PKG_MANAGER" == "brew" ]]; then
    echo -e "${BLUE}==> Using Apple Clang toolchain (macOS)${NC}"
    export CC=clang
    export CXX=clang++
    export FC=gfortran
    
    export LDFLAGS="-L$(brew --prefix libomp)/lib"
    export CPPFLAGS="-I$(brew --prefix libomp)/include"
    export CXXFLAGS="$CPPFLAGS"
else
    echo -e "${BLUE}==> Using GCC toolchain${NC}"
    export CC=gcc
    export CXX=g++
    export FC=gfortran
fi

# Fortran safety check
if ! command -v $FC &> /dev/null; then
    echo -e "${RED}[ERROR] Fortran compiler ($FC) not found. Required by PETSc.${NC}"
    exit 1
fi

# ------------------------------------------------------------
# PETSc & MPI DETECTION
# ------------------------------------------------------------
export PETSC_DIR="$INSTALL_DIR/petsc"

if ( [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" ] && \
     [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.dylib" ] && \
     [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a" ] ) || \
   [ "$FORCE_UPDATE" = true ]; then

    echo -e "${BLUE}==> Building PETSc ($PETSC_ARCH)...${NC}"
    [ ! -d "$PETSC_DIR" ] && git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"

    cd "$PETSC_DIR"
    DEBUG_FLAG=$([[ "$MODE_ARG" == "debug" ]] && echo "yes" || echo "no")
        
    # Intelligent MPI & Compiler handling
    if command -v mpicc &> /dev/null && command -v mpicxx &> /dev/null; then
        echo -e "${GREEN}    -> System MPI detected, skipping MPICH download.${NC}"
        MPI_OPTION="--with-mpi=1"
        
        if command -v mpif90 &> /dev/null; then
            MPI_FC="mpif90"
        elif command -v mpifort &> /dev/null; then
            MPI_FC="mpifort"
        else
            MPI_FC="$FC"
        fi
        COMPILER_OPTS="--with-cc=mpicc --with-cxx=mpicxx --with-fc=$MPI_FC"
    else
        echo -e "${YELLOW}    -> No System MPI detected, forcing MPICH download.${NC}"
        MPI_OPTION="--download-mpich"
        
        # On force GCC/Clang car PETSc va compiler MPICH lui-même
        COMPILER_OPTS="--with-cc=$CC --with-cxx=$CXX --with-fc=$FC"
    fi

    ./configure $COMPILER_OPTS \
                --with-debugging=$DEBUG_FLAG \
                --download-fblaslapack \
                $MPI_OPTION \
                --with-shared-libraries=1 \
                --PETSC_ARCH="$PETSC_ARCH" || exit 1

    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" -j$(($NCORES > 1 ? $NCORES - 1 : 1)) || exit 1
    cd - >/dev/null
else
    echo -e "${GREEN}==> PETSc ($PETSC_ARCH) is already built.${NC}"
fi

# ------------------------------------------------------------
# SYMENGINE (detection with build mode check)
# ------------------------------------------------------------
echo -e "${BLUE}==> Locating SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
SYM_FOUND=false
SEARCH_DIRS=("$PARENT_DIR" "$INSTALL_DIR" "/usr/local")

# Note: `find` is kept for dev-env comfort, but limited in depth to avoid HPC lags.

if [ -n "${SymEngine_DIR:-}" ]; then
    if [ -f "$SymEngine_DIR/SymEngineConfig.cmake" ]; then

        CACHE_FILE="$SymEngine_DIR/CMakeCache.txt"
        if [ -f "$CACHE_FILE" ]; then
            DETECTED_MODE=$(grep -i "^CMAKE_BUILD_TYPE" "$CACHE_FILE" | cut -d "=" -f 2 | tr '[:upper:]' '[:lower:]' | tr -d ' ')
            if [ -z "$DETECTED_MODE" ]; then
                echo -e "${YELLOW}[WARN] Unable to detect SymEngine build type (multi-config?)${NC}"
            elif [ "$DETECTED_MODE" != "$MODE_ARG" ]; then
                echo -e "${RED}[ERROR] SymEngine build type mismatch (expected $MODE_ARG, got $DETECTED_MODE)${NC}"
                exit 1
            fi
        fi

        echo -e "${GREEN}[OK] Using user-provided SymEngine: $SymEngine_DIR${NC}"
        SYM_FOUND=true
    else
        echo -e "${RED}[ERROR] Provided SymEngine_DIR is invalid: $SymEngine_DIR${NC}"
        exit 1
    fi
fi

if [ "$SYM_FOUND" = false ]; then
    while IFS= read -r FOUND_SYM; do
        SYM_CONFIG_DIR=$(dirname "$FOUND_SYM")
        CACHE_FILE="$SYM_CONFIG_DIR/CMakeCache.txt"

        if [ -f "$CACHE_FILE" ]; then
            DETECTED_MODE=$(grep -i "^CMAKE_BUILD_TYPE" "$CACHE_FILE" | cut -d "=" -f 2 | tr '[:upper:]' '[:lower:]' | tr -d ' ')
            [ "$DETECTED_MODE" != "$MODE_ARG" ] && continue
        else
            if [ ! -f "$SYM_CONFIG_DIR/SymEngineTargets-${MODE_ARG}.cmake" ] && \
            !( [ "$MODE_ARG" == "release" ] && [ -f "$SYM_CONFIG_DIR/SymEngineTargets-noconfig.cmake" ] ); then
                continue
            fi
        fi

        export SymEngine_DIR="$SYM_CONFIG_DIR"
        echo -e "${GREEN}[OK] SymEngine detected at: $SymEngine_DIR${NC}"
        SYM_FOUND=true
        break
    done < <(find "${SEARCH_DIRS[@]}" -maxdepth 4 -name "SymEngineConfig.cmake" 2>/dev/null)
fi

if [ "$SYM_FOUND" = false ] || [ "$FORCE_UPDATE" = true ]; then
    echo -e "${BLUE}==> Building SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
    SYM_SRC_DIR="$INSTALL_DIR/symengine"
    [ ! -d "$SYM_SRC_DIR" ] && git clone https://github.com/symengine/symengine "$SYM_SRC_DIR"

    export SymEngine_DIR="$SYM_SRC_DIR/build_${MODE_ARG}"
    mkdir -p "$SymEngine_DIR" && cd "$SymEngine_DIR"

    # Forcing toolchain consistency
    cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE \
             -DCMAKE_C_COMPILER="$CC" \
             -DCMAKE_CXX_COMPILER="$CXX" \
             -DBUILD_SHARED_LIBS=ON \
             -DWITH_SYMENGINE_RCP=OFF \
             -DWITH_LLVM=ON \
             -Dln_depend_mode=OFF || {
                echo -e "${RED}[ERROR] SymEngine config failed.${NC}"
                exit 1
             }

    make -j$(($NCORES > 1 ? $NCORES - 1 : 1)) || exit 1
    cd - >/dev/null
fi

# ------------------------------------------------------------
# CMAKE USER PRESETS
# ------------------------------------------------------------
echo -e "${BLUE}==> Finalizing CMake configuration${NC}"
USER_PRESET_FILE="$PROJECT_ROOT/CMakeUserPresets.json"

if [[ "$SymEngine_DIR" == *"$INSTALL_DIR/symengine"* ]]; then
    DYN_SYM_DIR="$INSTALL_DIR/symengine/build_\${presetName}"
else
    DYN_SYM_DIR="$SymEngine_DIR"
fi

if [ -f "$USER_PRESET_FILE" ] && [ "$FORCE_UPDATE" = false ]; then
    echo -e "${YELLOW}[WARNING] $USER_PRESET_FILE already exists.${NC}"
    echo -e "${BLUE}Please ensure your local preset paths match:${NC}"
    echo -e "  - PETSC_DIR:     ${GREEN}$PETSC_DIR${NC}"
    echo -e "  - SymEngine_DIR: ${GREEN}$DYN_SYM_DIR${NC}"
else
    cat <<EOF > "$USER_PRESET_FILE"
{
  "version": 3,
  "configurePresets": [
    {
      "name": "user-paths",
      "hidden": true,
      "cacheVariables": {
        "PETSC_DIR": "$PETSC_DIR",
        "SymEngine_DIR": "$DYN_SYM_DIR"
      }
    },
    {
      "name": "debug",
      "inherits": ["shared-debug", "user-paths"],
      "displayName": "My Local Debug",
      "cacheVariables": {
        "CMAKE_C_COMPILER": "$PETSC_DIR/${OS_PREFIX}-c-debug/bin/mpicc",
        "CMAKE_CXX_COMPILER": "$PETSC_DIR/${OS_PREFIX}-c-debug/bin/mpicxx"
      }
    },
    {
      "name": "release",
      "inherits": ["shared-release", "user-paths"],
      "displayName": "My Local Release",
      "cacheVariables": {
        "CMAKE_C_COMPILER": "$PETSC_DIR/${OS_PREFIX}-c-opt/bin/mpicc",
        "CMAKE_CXX_COMPILER": "$PETSC_DIR/${OS_PREFIX}-c-opt/bin/mpicxx"
      }
    }
  ]
}
EOF
    echo -e "${GREEN}[OK] Presets generated with MPI compilers.${NC}"
fi

# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------
echo -e "${GREEN}====================================================${NC}"
echo -e "${GREEN} DEPENDENCIES & ENVIRONMENT READY (${CMAKE_BUILD_TYPE})${NC}"
echo -e "${GREEN}====================================================${NC}"
echo -e "${BLUE} cmake --preset $MODE_ARG${NC}"
echo -e "${GREEN}====================================================${NC}"