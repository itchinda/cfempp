#include <gtest/gtest.h>
#include <petsc.h>
#include "cfemutils.h"
#include <omp.h>

int main(int argc, char **argv) {
    // omp_set_num_threads(1);
    ::testing::InitGoogleTest(&argc, argv);
    
    cfem::Environment::initialize(argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    return result;
}