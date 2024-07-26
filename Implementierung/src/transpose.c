#include <stdlib.h>
#include <errno.h>
#include <stdio.h>

#include "radixsort.h"
#include "cs_matrix.h"
#include "transpose.h"

/**
 * Computes the column pointers using the row pointers of a CSC matrix
 * @param indices       sorted row pointer array
 * @param dest          column pointer array
 * @param numCols       Number of columns in the matrix
 * @param valueCount    Number of non-zero values in the matrix
 * @return
 */
int computeColPtr(uint64_t* indices, uint64_t* dest, uint64_t numCols, uint64_t valueCount) {
    for (uint64_t i = 0; i <= numCols; ++i) {
        dest[i] = 0;
    }
    for (uint64_t i = 0; i < valueCount; i++) {
        if (indices[i]  < numCols) {
            dest[indices[i] + 1]++;
        }
    }
    for (uint64_t i = 0; i < numCols; i++) {
        dest[i + 1] += dest[i];
    }
    return 0;
}

/**
 * Transposes a given CSC matrix
 * @param a     matrix to transpose
 * @param a_t   transposed matrix
 * @return
 */
int transpose(struct cscMatrixTranspose* a, struct cscMatrix* a_t) {
    a_t->valueCount = a->valueCount;
    a_t->rows = a->columns;
    a_t->columns = a->rows;
    a_t->values = malloc((a_t->valueCount) * sizeof(float));
    a_t->rowIndices = malloc((a_t->valueCount) * sizeof(uint64_t));
    a_t->colPtr = malloc((a_t->columns+1) * sizeof(uint64_t));

    if (!a_t->values || !a_t->rowIndices || !a_t->colPtr) {
        free(a_t->values);
        free(a_t->rowIndices);
        free(a_t->colPtr);
        errno = ENOMEM;
        return 0;
    }
    for (uint64_t i = 0; i < a_t->valueCount; ++i) {
        a_t->values[i] = a->values[i];
        a_t->rowIndices[i] = a->colIndices[i];
    }

    errno = 0;
    radixSort(a->rowIndices, a_t->values, a_t->rowIndices, a_t->valueCount);
    if (errno == ENOMEM) {
        return 0;
    }

    computeColPtr(a->rowIndices, a_t->colPtr, a_t->columns, a->valueCount);
    return 1;
}

