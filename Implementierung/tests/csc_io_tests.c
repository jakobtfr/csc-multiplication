#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "cs_matrix.h"
#include "csc_io.h"
#include "csc_io_tests.h"

void static create_test_file(const char* filename, const char* content) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Unable to create test file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "%s", content);
    fclose(file);
}


int test_getColIndices_example() {
    
    struct cscMatrixTranspose m = {4, 4, 5};
    uint64_t rI[] = {0,2,1,1,3};
    uint64_t cPtr[] = {0,2,3,4,5};
    float val[] = {5,0.5,6,1,3};
    m.rowIndices = malloc(5*sizeof(uint64_t));
    m.colPtr = malloc(5*sizeof(uint64_t));
    m.values = malloc(5*sizeof(float));
    if (!m.rowIndices || !m.colPtr || !m.values) {
        free(m.rowIndices);
        free(m.colPtr);
        free(m.values);
        return 0;
    }

    for (int i = 0; i < 5; ++i) {
        m.rowIndices[i] = rI[i];
        m.colPtr[i] = cPtr[i];
        m.values[i] = val[i];
    }

    uint64_t cI[5];
    get_col_indices(&m, cI);
    uint64_t expected[] = {0,0,1,2,3};
    for (int i = 0; i < 5; ++i) {
        if (cI[i] != expected[i]) {
            printf("test_getColIndices_example failed.\nExpected: ");
            printUint64Vector(expected, 5);
            printf("But was: ");
            printUint64Vector(cI, 5);
            return 0;
        }
    }

    free(m.rowIndices);
    free(m.colPtr);
    free(m.values);
    return 1;
}
    

int test_parseCSCFile() {
    create_test_file("testMatrixA.txt", "4,4\n5,0.5,6,1,3\n0,2,1,1,3\n0,2,3,4,5\n");
    create_test_file("testMatrixB.txt", "4,4\n5,0.5,6,1,3\n0,2,1,1,3\n0,2,3,4,5\n");

    struct cscMatrixTranspose* matrixA = malloc(sizeof(struct cscMatrixTranspose));
    struct cscMatrix* matrixB = malloc(sizeof(struct cscMatrix));

    parse_csc_file("testMatrixA.txt", "testMatrixB.txt", matrixA, matrixB);

    float expectedA_values[] = {5,0.5,6,1,3};
    uint64_t expectedA_rowIndices[] = {0,2,1,1,3};
    uint64_t expectedA_colPtr[] = {0,2,3,4,5};
    
    float expectedB_values[] = {5,0.5,6,1,3};
    uint64_t expectedB_rowIndices[] = {0,2,1,1,3};
    uint64_t expectedB_colPtr[] = {0,2,3,4,5};

    if (matrixA->rows != 4 || matrixA->columns != 4 || matrixB->rows != 4 || matrixB->columns != 4) {
        printf("Test failed: Dimension mismatch\n");
        printf("Dimension mismatch: cannot multiply %ld by %lu matrix by a %ld by %ld matrix.\n",
                    matrixA->rows, matrixA->columns, matrixB->rows, matrixB->columns);
        return 0;
    }

    if (!cmp_float_vec_eq(matrixA->values, expectedA_values, matrixA->valueCount) ||
        !cmp_uint64_t_vec_eq(matrixA->rowIndices, expectedA_rowIndices, 5) ||
        !cmp_uint64_t_vec_eq(matrixA->colPtr, expectedA_colPtr, 5)) {
        printFloatVector(matrixA->values, matrixA->valueCount);
        printUint64Vector(matrixA->rowIndices, 5);
        printUint64Vector(matrixA->colPtr, 5);
        printf("Test failed: matrixA content mismatch\n");
        return 0;
    }

    if (!cmp_float_vec_eq(matrixB->values, expectedB_values, matrixB->valueCount) ||
        !cmp_uint64_t_vec_eq(matrixB->rowIndices, expectedB_rowIndices, 5) ||
        !cmp_uint64_t_vec_eq(matrixB->colPtr, expectedB_colPtr, 5)) {
        printFloatVector(matrixB->values, matrixB->valueCount);
        printUint64Vector(matrixB->rowIndices, 5);
        printUint64Vector(matrixB->colPtr, 5);
        printf("Test failed: matrixB content mismatch\n");
        return 0;
    }

    printf("\nTest passed, matrices are correctly parsed.\n");
    printf("%"PRIu64",%"PRIu64"\n", matrixA->rows, matrixA->columns);
    printFloatVector(matrixA->values, matrixA->valueCount);
    printUint64Vector(matrixA->rowIndices, 5);
    printUint64Vector(matrixA->colPtr, 5);
    printf("%"PRIu64",%"PRIu64"\n", matrixB->rows, matrixB->columns);
    printFloatVector(matrixB->values, matrixB->valueCount);
    printUint64Vector(matrixB->rowIndices, 5);
    printUint64Vector(matrixB->colPtr, 5);

    free_transpose_matrix(matrixA);
    free_csc_matrix(matrixB);
    remove("testMatrixA.txt");
    remove("testMatrixB.txt");

    return 1;
}


int test_result_to_file_fixed() {
    struct cscMatrix expected = {4, 4, 5, NULL, NULL, NULL};
    uint64_t rI[] = {0,2,1,1,3};
    uint64_t cPtr[] = {0,2,3,4,5};
    float val[] = {5,0.5,6,1,3};
    expected.rowIndices = malloc(5*sizeof(uint64_t));
    expected.colPtr = malloc(5*sizeof(uint64_t));
    expected.values = malloc(5*sizeof(float));
    if (!expected.rowIndices || !expected.colPtr || !expected.values) {
        free(expected.rowIndices);
        free(expected.colPtr);
        free(expected.values);
        return 0;
    }

    for (int i = 0; i < 5; ++i) {
        expected.rowIndices[i] = rI[i];
        expected.colPtr[i] = cPtr[i];
        expected.values[i] = val[i];
    }
    
    if (!expected.values || !expected.rowIndices || !expected.colPtr) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    result_to_file(&expected, "output1.txt");

    struct cscMatrix result;
    struct cscMatrixTranspose foo;

    parse_csc_file("output1.txt", "output1.txt", &foo, &result);

    if (5 != result.valueCount || 4 != result.rows || 4 != result.columns ||
        !cmp_float_vec_eq(result.values, val, 5) || !cmp_uint64_t_vec_eq(result.rowIndices, rI, 5)
        || !cmp_uint64_t_vec_eq(result.colPtr, cPtr, 5))
    {
        printf("Results of test_result_to_file_fixed:\n");
        printf("Test: testResultToFile failed.\n");
        printf("%"PRIu64",%"PRIu64"\n", result.rows, result.columns);
        printFloatVector(result.values, result.valueCount);
        printUint64Vector(result.rowIndices, 5);
        printUint64Vector(result.colPtr, 5);
        return 0;
    }

    printf("Results of test_result_to_file_fixed:");
    printf("\nTest passed: matrices are correctly written to the file.\n");
    printf("%"PRIu64",%"PRIu64"\n", result.rows, result.columns);
    printFloatVector(result.values, result.valueCount);
    printUint64Vector(result.rowIndices, 5);
    printUint64Vector(result.colPtr, 5);

    free(result.colPtr);
    free(result.rowIndices);
    free(result.values);
    free(foo.colIndices);
    free(foo.colPtr);
    free(foo.rowIndices);
    free(foo.values);
    free(expected.colPtr);
    free(expected.rowIndices);
    free(expected.values);

    remove("output1.txt");
    return 1;
}

int test_result_to_file_random(uint64_t maxSize, uint64_t minSize) {
    struct cscMatrix matrix = {0};
    matrix.rows = minSize + rand() % (maxSize - minSize);    
    matrix.columns = minSize + rand() % (maxSize - minSize);

    generate_csc_matr_rand(&matrix, 10, 3);

    if (!testVectorPurity(matrix.values, matrix.valueCount)) printf("matrix contains zeros");

    result_to_file(&matrix, "output2.txt");

    FILE* file = fopen("output2.txt", "r");
    if (!file) {
        perror("Unable to open file for reading");
        return 0;
    }

    struct cscMatrix result;
    struct cscMatrixTranspose foo;

    parse_csc_file("output2.txt", "output2.txt", &foo, &result);

    if (!cmp_csc_eq(&result, &matrix))
    {
        printf("\nResults of test_result_to_file_random:\n");
        printf("Test failed: matrices are not correctly written to the file.\n");
        return 0;
    }

    printf("\nResults of test_result_to_file_random:\n");
    printf("Test passed: matrices are correctly written to the file.\n");

    free(result.colPtr);
    free(result.rowIndices);
    free(result.values);
    free(foo.colIndices);
    free(foo.colPtr);
    free(foo.rowIndices);
    free(foo.values);
    free(matrix.colPtr);
    free(matrix.rowIndices);
    free(matrix.values);
    
    remove("output2.txt");
    return 1;
}

