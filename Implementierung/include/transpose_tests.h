#ifndef TRANSPOSE_TESTS_H
#define TRANSPOSE_TESTS_H

#include <stdint.h>

int test_computeColPtr_fixed();

int test_radixsort_fixed();

int test_radixsort();

int test_transpose_fixed();

int test_transpose_rand_fixed();

void test_transpose_rand_fixed_large();

int test_transpose_performance_fixed();

int test_transpose_rand(uint64_t minSize, uint64_t maxSize);

#endif
