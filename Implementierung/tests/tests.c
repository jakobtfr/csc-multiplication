#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include "cs_matrix.h"
#include "matrix_mul_tests.h"
#include "csc_io_tests.h"
#include "transpose_tests.h"
#include "matrix_mul.h"

static void print_runtime(clock_t start, clock_t end) {
    double seconds = ((double) end - start) / CLOCKS_PER_SEC;
    if (seconds > 60) {
        unsigned minutes = ((unsigned) seconds)/60;
        unsigned remSec = ((unsigned) seconds) % 60;
        printf("Runtime: %u min %u s\n", minutes, remSec);
    } else {
        printf("Runtime: %.1f s\n", seconds);
    }
}

static void test_randomness(unsigned k, unsigned n, unsigned c) {
    unsigned totals[k];
    memset(totals, 0, sizeof(totals));
    for (unsigned i = 0; i < n; ++i) {
        unsigned results[k];
        memset(results, 0, sizeof(results));
        for (unsigned j = 0; j < c; ++j) {
            results[rand() % k]++;
        }

        for(unsigned j = 0; j < k; ++j) {
            totals[j] += results[j];
        }
    }

    printf("\nRNG Test for integers between 0 and %u:\n\nAverage distribution"
            " for %u iterations of %u calls to rand() (value: count - "
            "percentage):\n", k, n, c);

    for(unsigned i = 0; i < k; ++i) {
        double avg = (double) totals[i]/n;
        printf("%d: %.2f - %.2f%%\n", i, avg, 100*(avg/c));
    }
}

static void test_csc_randomness(int n, int minD, int maxD) {
    struct cscMatrix matrices[n];

    printf("Testing average randomness of %d matrices with expected densities "
            "between %.0f%% and %.0f%%...\n", n, (((float) 100)/minD), 
            (((float) 100)/maxD));

    for(int i = 0; i < n; ++i) {
        uint64_t rows = 1 + rand() % 10;
        uint64_t columns = 1 + rand() % 10;
        matrices[i].rows = rows; 
        matrices[i].columns = columns;
        errno = 0;
        generate_csc_matr_rand(&matrices[i], minD, maxD);
        if (errno) {
            for (int j = 0; j < i; ++j) {
                free(matrices[j].colPtr);
                free(matrices[j].rowIndices);
                free(matrices[j].values);
            }
            perror("Memory error testing matrix randomness. Aborting.\n");
            return;
        }
    }

    float totalDensity = 0;
    float minDensity = 2;
    float maxDensity = -1;
    for (int i = 0; i < n; ++i) {
        int maxValCnt = matrices[i].rows * matrices[i].columns;
        float density = ((float) matrices[i].valueCount) / maxValCnt;
        totalDensity += density;
        if (density > maxDensity) maxDensity = density;
        if (density < minDensity) minDensity = density;
        free(matrices[i].colPtr);
        free(matrices[i].rowIndices);
        free(matrices[i].values);
    }
    printf("Average density: %.1f%%\nMaximum density: %.1f%%\n"
            "Minimum density: %.1f%%\n", 100*totalDensity/n, 100*maxDensity,
            100*minDensity);
}

int run_tests(){
    srand(time(NULL)); // Set rng seed to current time

    const int count = 22;
    int passed = 0;

    int res[count];

    // Do not change the order of the tests. The return value of some tests may 
    // have a different interpretation than the others, and they are identified
    // by their indices. If you want to add more tests, append them to the end
    // of the array (do not forget to update count) and handle its return value
    // at your discretion.

    res[0] = test_mul_id(matr_mult_csc);
    res[1] = test_getColIndices_example();
    res[2] = test_transpose_rand_fixed();                             // Index 2 is free to use
    res[3] = test_mul_example_squared(matr_mult_csc);
    res[4] = test_parseCSCFile();
    res[5] = test_result_to_file_fixed();
    res[6] = test_mul_rand_fixed(matr_mult_csc);
    res[7] = test_computeColPtr_fixed();
    res[8] = test_mul_rand_fixed_large(matr_mult_csc, 0);
    res[9] = test_mul_rand(matr_mult_csc, 10, 100, 1); // Console output of test can be modified
                                     // with the second argument. Check the docs
                                     // for more details
    res[10] = test_radixsort();
    res[11] = test_transpose_fixed();
    res[12] = test_result_to_file_random(300, 299);
    res[13] = test_mul_id(matr_mult_csc_V1);
    res[14] = test_mul_example_squared(matr_mult_csc_V1);
    res[15] = test_mul_rand_fixed(matr_mult_csc_V1);
    res[16] = test_mul_rand_fixed_large(matr_mult_csc_V1,0);
    res[17] = test_mul_rand(matr_mult_csc_V1, 10, 100, 1);
    res[18] = 0;
    res[19] = test_random_matrix_generation(10, 300);
    res[20] = test_mul_v2_id();
    res[21] = test_mul_rand_v2(10, 100, 1);

    for (int i = 0; i < count; ++i) {
        passed += (res[i] == 1);
    }

    printf("\nTest Results\n\n");
    printf("Summary:\n");
    for (int i = 0; i < count; ++i) {
        switch (i) {
            case 9:
            case 17:
            case 22:
                printf("Test %d: %s\n", i, res[i] ? "successful execution"
                        : "memory error encountered");
                break;
            case 20:
                if (res[i] == 1) printf("Test 20: passed\n");
                else printf("Test 20: failed with error code %d\n", res[i]);
                break;
            default:
                printf("Test %d: %s\n", i, (res[i] ? "passed" : "failed"));
                break;
        }
    }
    printf("\nPassed %d of %d tests.\n", passed, count);
    return passed == count;

}
// profiling sandbox. Use it to test performance however you want (i.e., 
// code in this function is not meant to be saved long-term)
int perf() { 
    logData = 1;
    return test_mul_rand(matr_mult_csc_V1, 1900, 2000, 1);
} 

int main () {
    int res;
    srand(time(NULL)); // Seed rng with current time - DO NOT COMMENT OUT
    // res = perf();
    res = run_tests();

    return 0;
}
