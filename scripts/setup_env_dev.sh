#!/bin/bash
# ============================================================
# CFEM++ INTERNAL BOOTSTRAP
# - Auto-installs dependencies & Builds PETSc/SymEngine
# ============================================================

GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo -e "${RED}[ERROR] Run with 'source scripts/setup_internal.sh [release|debug]'${NC}"
    exit 1
fi

echo -e "${BLUE}==> CFEM++ INTERNAL BOOTSTRAP${NC}"

# --- MODE ---
MODE_ARG=$(echo "${1:-release}" | tr '[:upper:]' '[:lower:]')
CMAKE_BUILD_TYPE=$([[ "$MODE_ARG" == "debug" ]] && echo "Debug" || echo "Release")
PETSC_ARCH=$([[ "$MODE_ARG" == "debug" ]] && echo "arch-linux-c-debug" || echo "arch-linux-c-opt")

# --- PATHS ---
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARENT_DIR="$(cd "$PROJECT_ROOT/.." && pwd)"
INSTALL_DIR="$HOME/software"
mkdir -p "$INSTALL_DIR"

# --- SYSTEM DEPS ---
echo -e "${BLUE}==> Verifying system dependencies${NC}"
sudo apt-get update >/dev/null 2>&1
sudo apt-get install -y build-essential gfortran cmake git pkg-config \
    libgmp-dev libmpfr-dev libmpc-dev libzstd-dev \
    libcurl4-openssl-dev llvm-dev libomp-dev >/dev/null 2>&1

# --- PETSc ---
export PETSC_DIR="$INSTALL_DIR/petsc"
if [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.so" ] && [ ! -f "$PETSC_DIR/$PETSC_ARCH/lib/libpetsc.a" ]; then
    echo -e "${BLUE}==> Building PETSc ($PETSC_ARCH)...${NC}"
    [ ! -d "$PETSC_DIR" ] && git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
    
    cd "$PETSC_DIR"
    DEBUG_FLAG=$([[ "$MODE_ARG" == "debug" ]] && echo "yes" || echo "no")
    ./configure --with-debugging=$DEBUG_FLAG --download-fblaslapack --download-mpich --with-shared-libraries=1 --PETSC_ARCH="$PETSC_ARCH" || return 1
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" -j$(nproc) || return 1
    cd - >/dev/null
fi
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# --- SymEngine (Intelligent Cache Detection) ---
echo -e "${BLUE}==> Locating SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
SYM_FOUND=false
SEARCH_DIRS=("$PARENT_DIR" "$INSTALL_DIR" "/usr/local" "/usr/lib")

# Recherche globale des fichiers Config
ALL_FOUND_SYMS=$(find "${SEARCH_BASE_DIRS[@]}" -maxdepth 5 -name "SymEngineConfig.cmake" 2>/dev/null)

for FOUND_SYM in $ALL_FOUND_SYMS; do
    SYM_CONFIG_DIR=$(dirname "$FOUND_SYM")
    CACHE_FILE="$SYM_CONFIG_DIR/CMakeCache.txt"
    
    # Cas A : C'est un dossier de build local (CMakeCache.txt existe)
    if [ -f "$CACHE_FILE" ]; then
        DETECTED_MODE=$(grep -i "^CMAKE_BUILD_TYPE" "$CACHE_FILE" | cut -d "=" -f 2 | tr '[:upper:]' '[:lower:]' | tr -d ' ')
        if [ "$DETECTED_MODE" == "$MODE_ARG" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] SymEngine ($DETECTED_MODE) detected in build tree at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
        
    # Cas B : C'est une installation système globale (make install)
    else
        # On vérifie si l'installation contient explicitement la cible pour notre mode
        if [ -f "$SYM_CONFIG_DIR/SymEngineTargets-${MODE_ARG}.cmake" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] System SymEngine (matches $MODE_ARG) detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
            
        # Parfois, CMake installe sans config spécifique ("noconfig"), souvent assimilé à un Release
        elif [ "$MODE_ARG" == "release" ] && [ -f "$SYM_CONFIG_DIR/SymEngineTargets-noconfig.cmake" ]; then
            export SymEngine_DIR="$SYM_CONFIG_DIR"
            echo -e "${GREEN}[OK] System SymEngine (noconfig -> Release) detected at: $SymEngine_DIR${NC}"
            SYM_FOUND=true; break
        fi
    fi
done

if [ "$SYM_FOUND" = false ]; then
    echo -e "${BLUE}==> Building SymEngine ($CMAKE_BUILD_TYPE)...${NC}"
    SYM_SRC_DIR="$PARENT_DIR/symengine"
    [ ! -d "$SYM_SRC_DIR" ] && git clone https://github.com/symengine/symengine "$SYM_SRC_DIR"
    
    export SymEngine_DIR="$SYM_SRC_DIR/build_${MODE_ARG}"
    mkdir -p "$SymEngine_DIR" && cd "$SymEngine_DIR"
    cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DWITH_LLVM=ON -Dln_depend_mode=OFF || return 1
    make -j$(nproc) || return 1
    cd - >/dev/null
fi

# --- BUILD CFEM ---
echo -e "${BLUE}==> Configuring CFEM++${NC}"
cd "$PROJECT_ROOT"
mkdir -p build_$MODE_ARG && cd build_$MODE_ARG
rm -f CMakeCache.txt

cmake .. -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -DPETSC_DIR="$PETSC_DIR" -DPETSC_ARCH="$PETSC_ARCH" -DSymEngine_DIR="$SymEngine_DIR" || return 1

echo -e "${GREEN}=============================================${NC}"
echo -e "${GREEN} CFEM++ INTERNAL READY (${CMAKE_BUILD_TYPE})${NC}"
echo -e "${GREEN} Run: cd build_$MODE_ARG && make -j${NC}"
echo -e "${GREEN}=============================================${NC}"