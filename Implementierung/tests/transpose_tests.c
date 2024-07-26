#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

#include "transpose_tests.h"
#include "transpose.h"
#include "cs_matrix.h"
#include "radixsort.h"
#include "csc_io.h"

/**
 * Prints transpose(input) = result in Wolfram language plain text
 *
 * @param result    Result of transpose(input)
 * @param input     Input matrix
 */
static void printExpr(struct cscMatrix *result, struct cscMatrix *input) {
    printf("transpose(");
    printCSCMatrix(input, 1);
    printf(") = ");
    printCSCMatrix(result, 1);
    printf("\n");
}

int test_radixsort_fixed() {
    struct cscMatrixTranspose A = {0};
    float valuesA[] = {5, 0.5f, 6, 1, 3};
    uint64_t rowIdxsA[] = {0, 2, 1, 1, 3};
    uint64_t colIdxsA[] = {0, 0, 1, 2, 3};
    uint64_t colPtrA[] = {0, 2, 3, 4, 5};

    uint64_t rowIdxsExp[] = {0, 1, 2, 0, 3};
    float valuesExp[] = {5, 6, 1, 0.5f, 3};
    uint64_t sortedRowIdx[] = {0, 1, 1, 2, 3};

    int res = 0;

    generateCSCMatrixTranspose(&A, 4, 4, 5, rowIdxsA, colIdxsA, colPtrA, valuesA);
    radixSort(rowIdxsA, valuesA, colIdxsA, 5);

    printf("\ntest_radixsort_fixed expected :\n Sorted Rows: \n");
    printUint64Vector(sortedRowIdx, 5);
    printf("Actual was:\n");
    printUint64Vector(rowIdxsA, 5);

    res = cmp_uint64_t_vec_eq(sortedRowIdx, rowIdxsA, 5);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Sorted Row mismatch. Test failed.\n");
    }

    printf("\ntest_radixsort_fixed expected :\n Values: \n");
    printFloatVector(valuesExp, 5);
    printf("Actual was:\n");
    printFloatVector(valuesA, 5);

    res = cmp_float_vec_eq(valuesA, valuesExp, 5);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Value mismatch. Test failed.\n");
    }

    printf("\ntest_radixsort_fixed expected :\n Row Indices: \n");
    printUint64Vector(rowIdxsExp, 5);
    printf("Actual was:\n");
    printUint64Vector(colIdxsA, 5);

    res = cmp_uint64_t_vec_eq(rowIdxsExp, colIdxsA, 5);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Column Index mismatch. Test failed.\n");
    }

    return res;
}

int test_radixsort() {
    float valuesA[] = {81.75, 29.3, 83.96, 19.21};
    uint64_t rowIdxsA[] = {2, 3, 0, 3};
    uint64_t colIdxsA[] = {0, 0, 2, 2};

    uint64_t colPtrA[] = {0, 2, 3, 4, 5};

    int rows = 4;
    int cols = 4;
    int valCount = 4;

    radixSort(rowIdxsA, valuesA, colIdxsA, valCount);
    printf("row indices:\n");
    printUint64Vector(rowIdxsA, valCount);
    printf("values:\n");
    printFloatVector(valuesA, 4);
    printf("colIdx:\n");
    printUint64Vector(colIdxsA, 4);
    return 1;
}

int test_computeColPtr_fixed() {
    uint64_t sortedRowIdx[] = {0, 1, 1, 2, 3};
    uint64_t result[5] = {0};
    uint64_t expected[] = {0, 1, 3, 4, 5};

    computeColPtr(sortedRowIdx, result, 4, 5);
    printf("\ntest_computeColPtr_fixed expected:\n");
    printUint64Vector(expected, 5);
    printf("Actual was:\n");
    printUint64Vector(result, 5);
    int res = cmp_uint64_t_vec_eq(expected, result, 5);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Column pointer mismatch. Test failed.\n");
    }
    return res;
}

int test_transpose_fixed() {
    struct cscMatrixTranspose matrA = {0};
    struct cscMatrix result = {0};
    struct cscMatrix expected = {0};

    float valuesExp[] = {5, 6, 1, 0.5f, 3};
    uint64_t rowIdxsExp[] = {0, 1, 2, 0, 3};
    uint64_t colPtrExp[] = {0, 1, 3, 4, 5};

    generate_csc_matr(&expected, 4, 4, 5, rowIdxsExp, colPtrExp, valuesExp);

    float valuesA[] = {5, 0.5f, 6, 1, 3};
    uint64_t rowIdxsA[] = {0, 2, 1, 1, 3};
    uint64_t colIdxsA[] = {0, 0, 1, 2, 3};
    uint64_t colPtrA[] = {0, 2, 3, 4, 5};

    generateCSCMatrixTranspose(&matrA, 4, 4, 5, rowIdxsA, colIdxsA, colPtrA, valuesA);

    transpose(&matrA, &result);
    printf("\ntest_tranpose_fixed expected:\n");
    printCSCMatrix(&expected, 0);
    printf("Actual was:\n");
    printCSCMatrix(&result, 0);
    int res = cmp_csc_eq(&expected, &result);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Transposed Matrix mismatch. Test failed.\n");
    }
    return res;
}

int test_transpose_rand_fixed() {
    struct cscMatrixTranspose a = {0};
    struct cscMatrix a_t = {0};
    struct cscMatrix exp = {0};
    struct cscMatrix cscMatrix = {0}; // input matrix as cscMatrix

    uint64_t rowIsA[] = {2, 3, 0, 3};
    uint64_t colPtrA[] = {0, 2, 2, 4};
    float aVals[] = {81.75, 29.3, 83.96, 19.21};
    uint64_t colIsA[] = {0, 0, 2, 2};

    uint64_t rowIsCSC[] = {2, 3, 0, 3};
    uint64_t colPtrCSC[] = {0, 2, 2, 4};
    float CSCVals[] = {81.75, 29.3, 83.96, 19.21};

    generate_csc_matr(&cscMatrix, 4, 3, 4, rowIsCSC, colPtrCSC, CSCVals);
    generateCSCMatrixTranspose(&a, 4, 3, 4, rowIsA, colIsA, colPtrA, aVals);

    float expVals[] = {83.96, 81.75, 29.3, 19.21};
    uint64_t rowIsExp[] = {2, 0, 0, 2};
    uint64_t colPtrExp[] = {0, 1, 1, 2, 4};

    generate_csc_matr(&exp, 3, 4, 4, rowIsExp, colPtrExp, expVals);
    printf("\ntest_transpose_rand expected:\n");
    printCSCMatrix(&exp, 0);
    transpose(&a, &a_t);
    printf("\nActual was:\n");
    printCSCMatrix(&a_t, 0);
    printf("\nWolfram Alpha comparison:\n");
    printExpr(&a_t, &cscMatrix);
    int res = cmp_csc_eq(&exp, &a_t);
    if (res) {
        printf("Test passed.\n");
    } else {
        printf("Transposed Matrix mismatch. Test failed.\n");
    }
    return res;
}

/**
 * This method is used for manual testing with wolfram alpha,
 * as the matrices are too large to test by hand
 */
void test_transpose_rand_fixed_large() {
    struct cscMatrixTranspose a = {0};
    struct cscMatrix cscMatrix = {0};
    struct cscMatrix a_t = {0};

    uint64_t rowIsA[] = {1, 1, 3, 4, 6, 8, 1, 4, 5, 9, 1, 5, 7, 2, 4, 8, 1, 2, 5, 1, 4, 6, 7, 8, 2, 8, 1,
                         4, 5, 6, 8, 0, 3, 4, 5, 6, 9, 1, 3, 6, 7, 8, 5, 7};
    uint64_t colPtrA[] = {0, 1, 6, 10, 13, 16, 19, 24, 26, 31, 37, 42, 44};
    float aVals[] = {47.4, 43.45, 7.13, 62.91, 68.25, 37.86, 11.59, 12.59,
                     85.33, 46.28, 84.38, 15.25, 34.79, 19.18, 30.52, 57.11, 1.82, 62.96,
                     42.79, 18.11, 67.43, 66.43, 79.42, 33.76, 55.35, 31.39, 33.64, 26.46,
                     68.14, 45.3, 69.69, 54.76, 74, 33.31, 49.28, 32.68, 81.29, 76.11,
                     61.45, 89.56, 89.2, 73.4, 87.92, 47.2};

    uint64_t *colIsA = malloc(sizeof(rowIsA) / sizeof(rowIsA[0]) * 
            sizeof(uint64_t));
    calculateColumnIndices(colPtrA, 12, colIsA);

    generateCSCMatrixTranspose(&a, 10, 12, 44, rowIsA, colIsA, colPtrA, aVals);

    transpose(&a, &a_t);

    uint64_t rowIsCSC[] = {1, 1, 3, 4, 6, 8, 1, 4, 5, 9, 1, 5, 7, 2, 4, 8, 1, 2, 5, 1, 4, 6, 7, 8, 2, 8, 1,
                           4, 5, 6, 8, 0, 3, 4, 5, 6, 9, 1, 3, 6, 7, 8, 5, 7};
    uint64_t colPtrCSC[] = {0, 1, 6, 10, 13, 16, 19, 24, 26, 31, 37, 42, 44};
    float CSCVals[] = {47.4, 43.45, 7.13, 62.91, 68.25, 37.86, 11.59, 12.59,
                       85.33, 46.28, 84.38, 15.25, 34.79, 19.18, 30.52, 57.11, 1.82, 62.96,
                       42.79, 18.11, 67.43, 66.43, 79.42, 33.76, 55.35, 31.39, 33.64, 26.46,
                       68.14, 45.3, 69.69, 54.76, 74, 33.31, 49.28, 32.68, 81.29, 76.11,
                       61.45, 89.56, 89.2, 73.4, 87.92, 47.2};


    generate_csc_matr(&cscMatrix, 10, 12, 44, rowIsCSC, colPtrCSC, CSCVals);

    printf("\ngiven: \n");
    printCSCMatrix(&cscMatrix, 0);

    printf("\nTransposed:\n");
    printCSCMatrix(&a_t, 0);

    printExpr(&a_t, &cscMatrix); // for manual testing with wolfram alpha

    free(colIsA);
}

int test_transpose_rand(uint64_t minSize, uint64_t maxSize) {
    struct cscMatrixTranspose a = {0};
    struct cscMatrix res = {0};
    uint64_t diff = maxSize - minSize;
    a.rows = minSize + rand() % diff;
    a.columns = minSize + rand() % diff;
    generateCSCMatrixTRand(&a, 10, 3);
    errno = 0;
    transpose(&a, &res);
    return !errno;
}
