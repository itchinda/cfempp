#!/bin/bash

# --- COLOR CONFIGURATION ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}==> [CFEM++] Initializing Environment Setup...${NC}"

# --- 1. PKG-CONFIG & SYSTEM DEPENDENCIES ---
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
                        libedit-dev libcurl4-openssl-dev llvm-dev libomp-dev

# --- 2. ENVIRONMENT VARIABLES PATHS ---
# Uses CI default paths if not already set in the shell
export PETSC_DIR=${PETSC_DIR:-$HOME/petsc}
export PETSC_ARCH=${PETSC_ARCH:-arch-linux-c-opt}
export SYMENGINE_INSTALL_DIR=${SYMENGINE_INSTALL_DIR:-$HOME/symengine_install}

echo -e "PETSC_DIR  : $PETSC_DIR"
echo -e "PETSC_ARCH : $PETSC_ARCH"

# --- 3. PETSC INSTALLATION (Automated if folder is missing) ---
if [ ! -d "$PETSC_DIR" ]; then
    echo -e "${BLUE}==> PETSc not detected. Cloning and building in $PETSC_DIR...${NC}"
    git clone -b release https://gitlab.com/petsc/petsc.git "$PETSC_DIR"
    cd "$PETSC_DIR"
    ./configure --with-debugging=no \
                --download-fblaslapack \
                --download-mpich \
                --with-shared-libraries=1 \
                --PETSC_ARCH="$PETSC_ARCH"
    make PETSC_DIR="$PETSC_DIR" PETSC_ARCH="$PETSC_ARCH" all
    cd -
else
    echo -e "${GREEN}[OK] PETSc found in $PETSC_DIR${NC}"
fi

# --- 4. SYMENGINE INSTALLATION (Automated if folder is missing) ---
if [ ! -d "$SYMENGINE_INSTALL_DIR" ]; then
    echo -e "${BLUE}==> SymEngine not detected. Cloning and building...${NC}"
    git clone https://github.com/symengine/symengine ~/symengine_source
    mkdir -p ~/symengine_source/build && cd ~/symengine_source/build
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX="$SYMENGINE_INSTALL_DIR" \
             -DWITH_SYMENGINE_RCP=ON \
             -DWITH_LLVM=ON \
             -Dln_depend_mode=OFF
    make -j$(nproc) install
    cd -
    rm -rf ~/symengine_source # Clean up source folder after installation
else
    echo -e "${GREEN}[OK] SymEngine found in $SYMENGINE_INSTALL_DIR${NC}"
fi

# --- 5. SESSION ENVIRONMENT UPDATE ---
export PATH="$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

# --- 6. CFEM++ PROJECT CONFIGURATION ---
echo -e "${BLUE}==> Configuring CFEM++ Project with CMake...${NC}"
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DPETSC_DIR="$PETSC_DIR" \
         -DPETSC_ARCH="$PETSC_ARCH" \
         -DSymEngine_DIR="$SYMENGINE_INSTALL_DIR/lib/cmake/symengine"

echo -e "${GREEN}============================================="
echo -e "  CFEM++ Environment is READY!"
echo -e "  Run 'cd build && make -j' to compile."
echo -e "=============================================${NC}"
