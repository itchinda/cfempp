//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#include <gtest/gtest.h>
#include <petsc.h>
#include "cfemutils.h"
#include <omp.h>

int main(int argc, char **argv) {
    omp_set_num_threads(1);
    ::testing::InitGoogleTest(&argc, argv);
    
    cfem::Environment::initialize(argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    return result;
}