#ifndef MAT_MUL_TESTS_H
#define MAT_MUL_TESTS_H

#include <stdint.h>

int test_mul_id(void (*mul_fun)(const void*, const void*, void*));

int test_mul_invalid_dims(void (*mul_fun)(const void*, const void*, void*));

int test_mul_example_squared(void (*mul_fun)(const void*, const void*, void*));

int test_mul_rand_fixed(void (*mul_fun)(const void*, const void*, void*));

int test_mul_rand_fixed_large(void (*mul_fun)(const void*, const void*, void*),
        unsigned verbose);

/**
 * Generates two random (in general non-square) matrices and computes their 
 * product.
 *
 * @param minSize       The minimum amount of rows and columns of the inputs
 * @param maxSize       The maximum amount of rows and columns of the inputs
 * @param verbosity     How much data the test outputs.
 *                      Verbosity 0:    Single-line success message
 *                      Verbosity 1:    Dimensions and density of all three 
 *                                      matrices, as well as runtime
 *                      Verbosity 2:    Verbosity 1 + whitespace+newline 
 *                                      representation of each matrix
 *                      Verbosity 3:    Verbosity 2 + the Wolfram language 
 *                                      representation of each matrix
 *                      Verbosity 4:    Verbosity 3 + the Wolfram language
 *                                      representation of the expression
 *                                      transpose(A) * B = Result
 *                      Verbosity 5:    Verbosity 4 + all the metadata
 *                                      of each matrix
 *
 * @return 1 if there were no memory errors during execution, 0 otherwise
 * 
 * It must hold: 0 < minSize <= maxSize. Otherwise, errno will be set to EINVAL
 *
 */
int test_mul_rand(void (*mul_fun)(const void*, const void*, void*)
        ,uint64_t minSize, uint64_t maxSize, int verbosity);

int test_random_matrix_generation(uint64_t minSize, uint64_t maxSize);

int test_mul_v2_id();

int test_mul_rand_v2(uint64_t minSize, uint64_t maxSize, int verbosity);

#endif
