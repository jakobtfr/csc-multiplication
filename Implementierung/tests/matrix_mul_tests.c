#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cs_matrix.h"
#include "matrix_mul.h"
#include "matrix_mul_tests.h"


static void printTestResult(const char* testName, struct cscMatrix* expected, 
        struct cscMatrix* actual) {
    printf("\nResults of test %s:\n", testName);
    printf("Test expected:\n");
    printCSCMatrix(expected, 0);
    printf("Result was:\n");
    printCSCMatrix(actual, 0);
}

static void printRandomTestResult(struct cscMatrix* a, struct cscMatrix* b,
        struct cscMatrix* result, unsigned verbosity, uint64_t maxSize, 
        int transpose) {
    if (verbosity > 5) {
        errno = EINVAL;
        perror("Random matrix-matrix product: verbosity must be an integer "
                "between 0 and 5. See the documentation for more information.");
        return;
    }
    printf("\nResult of random matrix-matrix multiplication A*B with maximum "
            "size %lu:\n", maxSize);
    if (verbosity == 0) {
        printf("Operation successful.\n");
        return;
    }
    if (verbosity == 1) {
        printf("Matrix A\n");
        printf("Dimensions: %lu rows and %lu columns (internal representation)\n"
                , a->rows, a->columns);
        printf("Density: %.2f%%\n\n", 
                ((float) a->valueCount*100)/(a->rows*a->columns));
        printf("Matrix B\n");
        printf("Dimensions: %lu rows and %lu columns (actual)\n"
                , b->rows, b->columns);
        printf("Density: %.2f%%\n\n", 
                ((float) b->valueCount*100)/(b->rows*b->columns));
        printf("Result\n");
        printf("Dimensions: %lu rows and %lu columns (actual)\n"
                , result->rows, result->columns);
        printf("Density: %.2f%%\n", 
                ((float) result->valueCount*100)/(result->rows*b->columns));
        return;
    }

    printf("Matrix A was:\n\n");
    printCSCMatrix(a, 0);
    printf("\nDimensions: %lu rows and %lu columns (internal representation)\n"
            , a->rows, a->columns);
    printf("Density: %.2f%%\n\n", ((float) a->valueCount*100)/(a->rows*a->columns));

    printf("Matrix B was:\n\n");
    printCSCMatrix(b, 0);
    printf("\nDimensions: %lu rows and %lu columns (actual)\n"
            , b->rows, b->columns);
    printf("Density: %.2f%%\n\n", ((float) b->valueCount*100)/(b->rows*b->columns));
    
    printf("Result was:\n\n");
    printCSCMatrix(result, 0);
    printf("\nDimensions %lu rows and %lu columns (actual)\n"
            , result->rows, result->columns);
    printf("Density: %.2f%%\n\n", 
            ((float) result->valueCount*100)/(result->rows*b->columns));
    if (verbosity == 2) return;

    printf("Wolfram language representation of matrices\n\n");
    printf("Matrix A\n");
    printCSCMatrix(a, 1);

    printf("\n\nMatrix B\n");
    printCSCMatrix(b, 1);

    printf("\n\nResult\n");
    printCSCMatrix(result, 1);

    if (verbosity == 3) return;

    
    printf("\n\nA * B - Result\n");
    if (transpose) printf("transpose(");
    printCSCMatrix(a, 1);
    if (transpose) printf(")");
    printf("*");
    printCSCMatrix(b, 1);
    printf("-");
    printCSCMatrix(result, 1);
    printf("\n");

    if (verbosity == 4) return;

    printf("Metadata of A:\n");
    printCSCMatrixMetadata(a);
    printf("\nMetadata of B:\n");
    printCSCMatrixMetadata(b);
    printf("\nMetadata of Result:\n");
    printCSCMatrixMetadata(result);
}

static int compareResultExpected(struct cscMatrix* res, struct cscMatrix* ex) {
    int cmpResult = cmp_csc_eq(res, ex);
    if (!cmpResult) {
        printf("cmp_csc_eq returned false. Test failed.\n");
    } else {
        printf("cmp_csc_eq returned true; matrices are within tolerable"
                " epsilon. Test passed.\n");
    }
    return cmpResult;
}

int test_mul_id(void (*mul_fun)(const void*, const void*, void*)) {
    struct cscMatrix matrA = {0};
    struct cscMatrix matrB = {0};
    struct cscMatrix result = {0};
    struct cscMatrix expected = {0};

    float valuesA[] = {5,0.5f,6,1,3};
    uint64_t rowIdxsA[] = {0,2,1,1,3};
    uint64_t colPtrA[] = {0,2,3,4,5};

    generate_csc_matr(&matrA, 4, 4, 5, rowIdxsA, colPtrA, valuesA);

    float valuesB[4];
    uint64_t rowIdxsB[4];
    uint64_t colPtrB[5];

    matrB.values = valuesB;
    matrB.rowIndices = rowIdxsB;
    matrB.colPtr = colPtrB;
    generate_csc_id(&matrB, 4);

    errno = 0;
    (*mul_fun)(&matrA, &matrB, &result);
    if (errno != 0) return 0;


    float exValues[] = {5,6,1,0.5,3};
    uint64_t exRowIdxs[] = {0,1,2,0,3};
    uint64_t exColPtr[] = {0,1,3,4,5};

    generate_csc_matr(&expected, 4, 4, 5, exRowIdxs, exColPtr, exValues);

    printTestResult("test_mul_id", &expected, &result);

    int testResult = compareResultExpected(&result, &expected);
    free(result.colPtr);
    free(result.values); 
    free(result.rowIndices);
    return testResult;
}

int test_mul_example_squared(void (*mul_fun)(const void*, const void*, void*)) {
    struct cscMatrix a = {0};
    struct cscMatrix b = {0};
    struct cscMatrix res = {0};
    struct cscMatrix x = {0};

    uint64_t rowIsA[] = {0, 1, 2, 0, 3};
    uint64_t colPtrA[] = {0, 1, 3, 4, 5 };
    float aVals[] = {5, 6, 1, 0.5, 3};

    generate_csc_matr(&a, 4, 4, 5, rowIsA, colPtrA, aVals);

    uint64_t rowIsB[] = {0, 2, 1, 1, 3};
    uint64_t colPtrB[] = {0, 2, 3, 4, 5};
    float bVals[] = {5, 0.5, 6, 1, 3};

    generate_csc_matr(&b, 4, 4, 5, rowIsB, colPtrB, bVals);

    uint64_t xRowIs[] = {0, 1, 2, 1, 1, 3};
    uint64_t xColPtr[] = {0, 3, 4, 5, 6};
    float xVals[] = {25, 0.5, 2.5, 36, 6, 9};

    generate_csc_matr(&x, 4, 4, 6, xRowIs, xColPtr, xVals);

    errno = 0;
    (*mul_fun)(&a, &b, &res);
    if (errno) return 0;

    printTestResult("test_mul_example_squared", &x, &res);
    
    int testResult = compareResultExpected(&res, &x);
    free(res.values);
    free(res.colPtr);
    free(res.rowIndices);
    return testResult;
}

int test_mul_rand_fixed(void (*mul_fun)(const void*, const void*, void*)) {
    struct cscMatrix a = {0};
    struct cscMatrix b = {0};
    struct cscMatrix res = {0};
    struct cscMatrix x = {0};

    uint64_t rowIsA[] = {2, 3, 0, 3};
    uint64_t colPtrA[] = {0, 2, 2, 4};
    float aVals[] = {81.75, 29.3, 83.96, 19.21};

    generate_csc_matr(&a, 4, 3, 4, rowIsA, colPtrA, aVals);

    uint64_t rowIsB[] = {0, 2, 2, 3, 0, 2, 3};
    uint64_t colPtrB[] = {0, 2, 3, 4, 5, 7};
    float bVals[] = {93.76, 19.22, 56.68, 35.56, 56.35, 96.7, 90.98};

    generate_csc_matr(&b, 4, 5, 7, rowIsB, colPtrB, bVals);

    uint64_t xRowIs[] = {0, 2, 0, 0, 2, 2, 0, 2};
    uint64_t xColPtr[] = {0, 2, 3, 5, 6, 8};
    float xVals[] = {1571.24, 7872.09, 4633.59, 1041.91, 683.108, 4731.15, 
        10570.9, 1747.73};

    generate_csc_matr(&x, 3, 5, 8, xRowIs, xColPtr, xVals);

    errno = 0;
    (*mul_fun)(&a, &b, &res);
    if (errno) {
        printf("test_mul_rand_fixed failed with errno value %d", errno);
        return 0;
    }

    printTestResult("test_mul_rand_fixed", &x, &res);
    int testResult = compareResultExpected(&res, &x);
    free(res.values);
    free(res.colPtr);
    free(res.rowIndices);
    return testResult;
}

/**
 * Tests product of (relatively) large statically-generated matrices.
 *
 * @param verbose   Set to 1 to display expected and result matrices, 
 *                  set to 0 to only display success/failure message
 * @return          1 if result equals expected, 0 otherwise
 */
int test_mul_rand_fixed_large(void (*mul_fun)(const void*, const void*, void*),
        unsigned verbose) {
    // Inputs and Expected in Wolfram language plain text:
    //
    // a: {{0,47.4,0,0,0,0,0,0,0,0},{0,43.45,0,7.13,62.91,0,68.25,0,37.86,0},{0,11.59,0,0,12.59,85.33,0,0,0,46.28},{0,84.38,0,0,0,15.25,0,34.79,0,0},{0,0,19.18,0,30.52,0,0,0,57.11,0},{0,1.82,62.96,0,0,42.79,0,0,0,0},{0,18.11,0,0,67.43,0,66.43,79.42,33.76,0},{0,0,55.35,0,0,0,0,0,31.39,0},{0,33.64,0,0,26.46,68.14,45.30,0,69.69,0},{54.76,0,0,74.0,33.31,49.28,32.68,0,0,81.29},{0,76.11,0,61.45,0,0,89.56,89.2,73.40,0},{0,0,0,0,0,87.92,0,47.20,0,0}} - 12*10, 36.7% Density
    //
    // b: {{0,0,0,76.90,62.59,0,0,0},{0,22.33,0,0,0,0,0,0},{0,0,34.59,0,0,0,0,0},{41.91,0,0,0,0,0,40.19,0},{0,0,0,0,0,0,26.63,0},{0,0,0,76.12,0,0,0,0},{4.3,58.93,0,0,0,0,0,66.99},{0,7.37,0,0,0,0,0,0},{0,0,0,0,0,39.9,91.14,0},{75.77,0,0,3.85,0,0,0,0}} - 10*8, 20% Density
    //
    // x: {{0, 1058.44, 0, 0, 0, 0, 0, 0}, {592.293, 4992.21, 0, 0, 0, 1510.61, 5412.41, 4572.07}, {3506.64, 258.805, 0, 6673.5, 0, 0, 335.272, 0}, {0, 2140.61, 0, 1160.83, 0, 0, 0, 0}, {0, 0, 663.436, 0, 0, 2278.69, 6017.75, 0}, {0, 40.6406, 2177.79, 3257.17, 0, 0, 0, 0}, {285.649, 4904.44, 0, 0, 0, 1347.02, 4872.55, 4450.15}, {0, 0, 1914.56, 0, 0, 1252.46, 2860.88, 0}, {194.79, 3420.71, 0, 5186.82, 0, 2780.63, 7056.18, 3034.65}, {9401.21, 1925.83, 0, 8275.2, 3427.43, 0, 3861.11, 2189.23}, {2960.48, 7634.71, 0, 0, 0, 2928.66, 9159.35, 5999.62}, {0, 347.864, 0, 6692.47, 0, 0, 0, 0}} - 12*8

    if (verbose > 1) {
        errno = EINVAL;
        perror("'verbose' can only have a value of 1 or 0");
        return 0;
    }
    struct cscMatrix a = {0};
    struct cscMatrix b = {0};
    struct cscMatrix res = {0};
    struct cscMatrix x = {0};

    uint64_t rowIsA[] = {1,1,3,4,6,8,1,4,5,9,1,5,7,2,4,8,1,2,5,1,4,6,7,8,2,8,1,
        4,5,6,8,0,3,4,5,6,9,1,3,6,7,8,5,7};
    uint64_t colPtrA[] = {0,1,6,10,13,16,19,24,26,31,37,42,44};
    float aVals[] = {47.4, 43.45, 7.13, 62.91, 68.25, 37.86, 11.59, 12.59, 
        85.33, 46.28, 84.38, 15.25, 34.79, 19.18, 30.52, 57.11, 1.82, 62.96, 
        42.79, 18.11, 67.43, 66.43, 79.42, 33.76, 55.35, 31.39, 33.64, 26.46,
        68.14, 45.3, 69.69, 54.76, 74, 33.31, 49.28, 32.68,  81.29, 76.11, 
        61.45, 89.56, 89.2, 73.4, 87.92, 47.2};

    generate_csc_matr(&a, 10, 12, 44, rowIsA, colPtrA, aVals);

    uint64_t rowIsB[] = {3, 6, 9, 1, 6, 7, 2, 0, 5, 9, 0, 8, 3, 4, 8, 6};
    uint64_t colPtrB[] = {0, 3, 6, 7, 10, 11, 12, 15, 16};
    float bVals[] = {41.91, 4.3, 75.77, 22.33, 58.93, 7.37, 34.59, 76.9, 76.12,
        3.85, 62.59, 39.9, 40.19, 26.63, 91.14, 66.99};

    generate_csc_matr(&b, 10, 8, 16, rowIsB, colPtrB, bVals);

    uint64_t xRowIs[] = {1,2,6,8,9,10,0,1,2,3,5,6,8,9,10,11,4,5,7,2,3,5,8,9,11,
        9,1,4,6,7,8,10,1,2,4,6,7,8,9,10,1,6,8,9,10};
    uint64_t xColPtr[] = {0, 6, 16, 19, 25, 26, 32, 40, 45};
    float xVals[] = {592.293, 3506.64, 285.649, 194.79, 9401.21, 2960.48, 
        1058.44, 4992.21, 258.805, 2140.61, 40.6406, 4904.44, 3420.71, 1925.83, 
        7634.71, 347.864, 663.436, 2177.79, 1914.56, 6673.5, 1160.83, 3257.17, 
        5186.82, 8275.2, 6692.47, 3427.43, 1510.61, 2278.69, 1347.02, 1252.46, 
        2780.63, 2928.66, 5412.41, 335.272, 6017.75, 4872.55, 2860.88, 7056.18,
        3861.11, 9159.35, 4572.07, 4450.15, 3034.65, 2189.23, 5999.62};

    generate_csc_matr(&x, 12, 8, 45, xRowIs, xColPtr, xVals);

    errno = 0;
    (*mul_fun)(&a, &b, &res);
    if (errno) {
        printf("test_mul_rand_fixed_large failed with errno value %d", errno);
        return 0;
    }

    if (verbose == 1) {
        printTestResult("test_mul_rand_fixed_large", &x, &res);
    } else {
        printf("\ntest_mul_rand_fixed_large: ");
    }
    
    int testResult = compareResultExpected(&res, &x);
    free(res.values);
    free(res.colPtr);
    free(res.rowIndices);
    return testResult;
}

int test_mul_rand(void (*mul_fun)(const void*, const void*, void*), 
        uint64_t minSize, uint64_t maxSize, int verbosity) {
    if (minSize == 0 || minSize > maxSize) {
        errno = EINVAL;
        perror("test_mul_rand: minSize must be greater than 0 and less or equal"
                " to maxSize");
        return 0;
    }

    struct cscMatrix a = {0};
    struct cscMatrix b = {0};
    struct cscMatrix res = {0};

    uint64_t diff = maxSize - minSize + 1;

    uint64_t aCols = minSize + rand() % diff;
    uint64_t abRows = minSize + rand() % diff;
    uint64_t bCols = minSize + rand() % diff;

    a.rows = b.rows = abRows;
    a.columns = aCols;
    b.columns = bCols;

    printf("\nBegin random matrix-matrix multiplication test\n");
    errno = 0;
    generate_csc_matr_rand(&a, 10, 3);
    if (errno) {
        perror("generate_csc_matr_rand A had a memory error.");
        return 0;
    }
    printf("Matrix A initialized successfully\n");

    generate_csc_matr_rand(&b, 10, 3);
    if (errno) {
        perror("generate_csc_matr_rand B had a memory error.");
        return 0;
    }
    printf("Matrix B initialized successfully\n");
    printf("Begin operation\n");
    clock_t start = clock();
    (*mul_fun)(&a, &b, &res);
    clock_t end = clock();
    if (errno) {
        perror("matr_mult_csc had a memory error.");
        return 0;
    }

    printRandomTestResult(&a, &b, &res, verbosity, maxSize, 1);
    double seconds = ((double) end - start) / CLOCKS_PER_SEC;
    if (seconds > 60) {
        unsigned minutes = ((unsigned) seconds)/60;
        unsigned remSec = ((unsigned) seconds) % 60;
        printf("Runtime: %u min %u s\n", minutes, remSec);
    } else {
        printf("Runtime: %.1f s\n", seconds);
    }

    int isPure = testVectorPurity(res.values, res.valueCount);
    int resVal = !errno && isPure;

    if (!isPure) {
        errno = EPERM;
        perror("Result of random matrix-matrix multiplication is not pure. "
                "Test failed");
    }

    free(a.values);
    free(a.rowIndices);
    free(a.colPtr);
    free(b.values);
    free(b.rowIndices);
    free(b.colPtr);
    free(res.values);
    free(res.rowIndices);
    free(res.colPtr);
    return resVal;
}

int test_random_matrix_generation(uint64_t minSize, uint64_t maxSize) {
    struct cscMatrix a = {0};
    uint64_t diff = maxSize - minSize;
    a.columns = minSize + rand() % diff;
    a.rows = minSize + rand() % diff;

    struct cscMatrixTranspose b = {0};
    b.columns = minSize + rand() % diff;
    b.rows = minSize + rand() % diff;

    generate_csc_matr_rand(&a, 10, 3);
    generateCSCMatrixTRand(&b, 10, 3);
    int resA = is_valid_csc(&a);
    if(resA != 1) {
        printf("test_random_matrix_generation failed: a is not valid\n");
        printCSCMatrixMetadata(&a);
        free(a.values);
        free(a.rowIndices);
        free(a.colPtr);
        free(b.values);
        free(b.rowIndices);
        free(b.colPtr);
        free(b.colIndices);
        return resA;
    }   
    int resB = is_valid_csc_transpose(&b);
    if(resB != 1) {
        printf("test_random_matrix_generation failed: b is not valid\n");
        printCSCMatrixMetadataT(&b);
        free(a.values);
        free(a.rowIndices);
        free(a.colPtr);
        free(b.values);
        free(b.rowIndices);
        free(b.colPtr);
        free(b.colIndices);
        return resB;
    }   
    free(a.values);
    free(a.rowIndices);
    free(a.colPtr);
    free(b.values);
    free(b.rowIndices);
    free(b.colPtr);
    free(b.colIndices);

    printf("test_random_matrix_generation passed: both matrices are valid\n");
    return 1;
}


int test_mul_v2_id() {
    struct cscMatrix matrA = {0};
    struct cscMatrix matrB = {0};
    struct cscMatrix result = {0};
    struct cscMatrix expected = {0};

    float valuesA[] = {5,0.5f,6,1,3};
    uint64_t rowIdxsA[] = {0,2,1,1,3};
    uint64_t colPtrA[] = {0,2,3,4,5};

    generate_csc_matr(&matrA, 4, 4, 5, rowIdxsA, colPtrA, valuesA);

    float valuesB[4];
    uint64_t rowIdxsB[4];
    uint64_t colPtrB[5];

    matrB.values = valuesB;
    matrB.rowIndices = rowIdxsB;
    matrB.colPtr = colPtrB;
    generate_csc_id(&matrB, 4);

    errno = 0;
    matr_mult_csc_V2(&matrA, &matrB, &result);
    if (errno != 0) return 0;


    float exValues[] = {5,0.5f,6,1,3};
    uint64_t exRowIdxs[] = {0,2,1,1,3};
    uint64_t exColPtr[] = {0,2,3,4,5};

    generate_csc_matr(&expected, 4, 4, 5, exRowIdxs, exColPtr, exValues);

    printTestResult("test_mul_id", &expected, &result);

    int testResult = compareResultExpected(&result, &expected);
    free(result.colPtr);
    free(result.values); 
    free(result.rowIndices);
    return testResult;
}


int test_mul_rand_v2(uint64_t minSize, uint64_t maxSize, int verbosity) {
    if (minSize == 0 || minSize > maxSize) {
        errno = EINVAL;
        perror("test_mul_rand: minSize must be greater than 0 and less or equal"
                " to maxSize");
        return 0;
    }

    struct cscMatrix a = {0};
    struct cscMatrix b = {0};
    struct cscMatrix res = {0};

    uint64_t diff = maxSize - minSize + 1;

    uint64_t aColsBRows = minSize + rand() % diff;
    uint64_t aRows = minSize + rand() % diff;
    uint64_t bCols = minSize + rand() % diff;

    a.columns = b.rows = aColsBRows;
    a.rows = aRows;
    b.columns = bCols;

    printf("\nBegin random matrix-matrix multiplication test\n");
    errno = 0;
    generate_csc_matr_rand(&a, 6, 2);
    if (errno) {
        perror("generate_csc_matr_rand A had a memory error.");
        return 0;
    }
    printf("Matrix A initialized successfully\n");

    generate_csc_matr_rand(&b, 6, 2);
    if (errno) {
        perror("generate_csc_matr_rand B had a memory error.");
        return 0;
    }
    printf("Matrix B initialized successfully\n");
    printf("Begin operation\n");
    clock_t start = clock();
    matr_mult_csc_V2(&a, &b, &res);
    clock_t end = clock();
    if (errno) {
        perror("matr_mult_csc had a memory error.");
        return 0;
    }

    printRandomTestResult(&a, &b, &res, verbosity, maxSize, 0);
    double seconds = ((double) end - start) / CLOCKS_PER_SEC;
    if (seconds > 60) {
        unsigned minutes = ((unsigned) seconds)/60;
        unsigned remSec = ((unsigned) seconds) % 60;
        printf("Runtime: %u min %u s\n", minutes, remSec);
    } else {
        printf("Runtime: %.1f s\n", seconds);
    }

    int isPure = testVectorPurity(res.values, res.valueCount);
    int resVal = !errno && isPure;

    if (!isPure) {
        errno = EPERM;
        perror("Result of random matrix-matrix multiplication is not pure. "
                "Test failed");
    }

    free(a.values);
    free(a.rowIndices);
    free(a.colPtr);
    free(b.values);
    free(b.rowIndices);
    free(b.colPtr);
    free(res.values);
    free(res.rowIndices);
    free(res.colPtr);
    return resVal;
}


