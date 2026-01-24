#include <gtest/gtest.h>
#include <petsc.h>

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    int result = RUN_ALL_TESTS();
    
    PetscFinalize();
    
    return result;
}