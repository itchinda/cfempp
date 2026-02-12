#include <gtest/gtest.h>
#include <petsc.h>
#include "cfemutils.h"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    cfem::Environment::initialize(argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    return result;
}