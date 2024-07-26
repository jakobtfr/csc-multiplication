#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <errno.h>

#include "cs_matrix.h"


int logData = 0;
                                                                                
void get_col_vector (float* restrict values, uint64_t* restrict indices, uint64_t len,
                     float* restrict vector, uint64_t* restrict subIndices, 
                     uint64_t start) {
    for (uint64_t i = 0; i < len; ++i) {
        vector[i] = values[start+i];
        subIndices[i] = indices[start+i];
    }
}
                                                                              
uint64_t extend_vector(float** values, uint64_t** indices, uint64_t current,
                       uint64_t max) {
    if (current == max) {
        errno = ENOMEM;
        return 0;
    }
    uint64_t doubleSize;
    uint64_t newSize;
    int of = __builtin_umull_overflow(2, current, 
            &doubleSize);
    newSize = (of || doubleSize > max) ? max : doubleSize;
    float* newVals = reallocarray(*values, newSize, sizeof(float));
    uint64_t* newInd = reallocarray(*indices, newSize, sizeof(uint64_t));
    if (!newVals || !newInd) {
        free(newVals);
        free(newInd);
        free(*values);
        free(*indices);
        errno = ENOMEM;
        return 0;
    }
    *values = newVals;
    *indices = newInd;
    return newSize;
}

struct cscRow get_row_vector(const struct cscMatrix* m, uint64_t index) {
    struct cscRow row = {
        .values = malloc(sizeof(float)),
        .colIndices = malloc(sizeof(uint64_t))
    };
    uint64_t size = 1;
    if (!row.values || !row.colIndices) {
        free(row.values);
        free(row.colIndices);
        row.values = 0;
        row.colIndices = 0;
        errno = ENOMEM;
        return row;
    }

    for (uint64_t i = 0; i < m->columns; ++i) {
        for (uint64_t j = m->colPtr[i]; j < m->colPtr[i + 1]; ++j) {
            if (m->rowIndices[j] != index) continue;

            if (row.valueCount >= size) {
                size = extend_vector(&row.values, &row.colIndices, row.valueCount,
                                     m->columns);
                if (!size) {
                    free(row.values);
                    free(row.colIndices);
                    row.values = 0;
                    row.colIndices = 0;
                    errno = ENOMEM;
                    return row;
                }
            }

            row.values[row.valueCount] = m->values[j];
            row.colIndices[row.valueCount++] = i;
            break;
        }
    }
    return row;
}

float scalar_prod(const float* restrict aVec, const float* restrict bVec, 
        const uint64_t* restrict aInd, const uint64_t* restrict bInd, 
        uint64_t aLen, uint64_t bLen) {

    float result = 0;

    for (uint64_t i = 0, j = 0; i < aLen && j < bLen;) {
        if (aInd[i] < bInd[j]) {
            while (i < aLen && aInd[i] < bInd[j]) ++i;

            if (i >= aLen) break;

        } else {
            while (j < bLen && bInd[j] < aInd[i]) ++j;

            if(j >= bLen) break;
        }

        if (aInd[i] == bInd[j]) {
            result += aVec[i++] * bVec[j++];
            }
    }

    return result;
}

float scalar_prod_in_place(const float* restrict aVec, 
        const float* restrict bVec, const uint64_t* restrict aInd, 
        const uint64_t* restrict bInd, uint64_t aStart, uint64_t aEnd, 
        uint64_t bStart, uint64_t bEnd) {
    float result = 0;

    for (uint64_t i = aStart, j = bStart; i < aEnd && j < bEnd;) {
        if (aInd[i] < bInd[j]) {
            while (i < aEnd && aInd[i] < bInd[j]) ++i;

            if (i >= aEnd) break;

        } else {
            while (j < bEnd && bInd[j] < aInd[i]) ++j;

            if(j >= bEnd) break;
        }

        if (aInd[i] == bInd[j]) {
            result += aVec[i++] * bVec[j++];
            }
    }

    return result;
}

int cmp_float_eq(float a, float b) {
    const float EPSILON = 1e-5; 

    if (a == b) return 1;
    float diff = a - b;
    if (diff < 0) diff = -diff;

    float absA = a > 0 ? a : -a, 
          absB = b > 0 ? b : -b;

    if (a == 0 || b == 0 || absA + absB < FLT_MIN) return diff < EPSILON;
    return diff/a < EPSILON;
}

int cmp_uint64_t_vec_eq(uint64_t* a, uint64_t* b, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}

int cmp_float_vec_eq(float* a, float* b, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) {
        if (!cmp_float_eq(a[i], b[i])) return 0;
    }
    return 1;
}

int cmp_csc_eq(struct cscMatrix* a, struct cscMatrix* b) {
    if (a->rows != b->rows) return 0;
    if (a->columns != b->columns) return 0;
    if (a->valueCount != b->valueCount) return 0;

    for (uint64_t i = 0; i < a->valueCount; ++i) {
        if (!cmp_float_eq(a->values[i], b->values[i])) return 0;
        if (a->rowIndices[i] != b->rowIndices[i]) return 0;
    }

    for(uint64_t i = 0; i < a->columns+1; ++i) {
        if (a->colPtr[i] != b->colPtr[i]) return 0;
    }

    return 1;
}

int generate_csc_id(struct cscMatrix* id, uint64_t n) {
    if (n < 1) return 0;
    id->columns = id->rows = n;
    id->colPtr[0] = 0;
    for (uint64_t i = 0; i < n; ++i) {
        id->values[i] = 1;
        id->colPtr[i+1] = i+1;
        id->rowIndices[i] = i;
    }
    return 1;
}

void generate_csc_matr(struct cscMatrix* dest, uint64_t rows, uint64_t columns,
                       uint64_t valueCount, uint64_t* rowIndices, uint64_t* colPtr,
                       float* values) {
    dest->rows = rows;
    dest->columns = columns;
    dest->valueCount = valueCount;
    dest->rowIndices = rowIndices;
    dest->colPtr = colPtr;
    dest->values = values;
}

void generateCSCMatrixTranspose(struct cscMatrixTranspose* dest, uint64_t rows, 
        uint64_t columns, uint64_t valueCount, uint64_t* rowIndices, 
        uint64_t* colIndices, uint64_t* colPtr, float* values) {
    dest->rows = rows;
    dest->columns = columns;
    dest->valueCount = valueCount;
    dest->rowIndices = rowIndices;
    dest->colIndices = colIndices;
    dest->colPtr = colPtr;
    dest->values = values;
}

void calculateColumnIndices(uint64_t *colPtr, uint64_t numCols, uint64_t *colInd) {
    for (uint64_t j = 0; j < numCols; j++) {
        for (uint64_t i = colPtr[j]; i < colPtr[j + 1]; i++) {
            colInd[i] = j;
        }
    }
}


void generate_csc_matr_rand(struct cscMatrix* dest, int minD_inv,
                            int maxD_inv) {

    dest->valueCount = 0;
    uint64_t maxVals;
    bool of = __builtin_umull_overflow(dest->rows, dest->columns, &maxVals);
    if (of) {
        errno = ERANGE;
        return;
    }
    float* vals = malloc(maxVals * sizeof(float));
    uint64_t* rowIndices = malloc(maxVals * sizeof(uint64_t));
    dest->colPtr = calloc(dest->columns+1, sizeof(uint64_t));
    if (!vals || !rowIndices || !dest->colPtr) {
        free(vals);
        free(rowIndices);
        free(dest->colPtr);
        errno = ENOMEM;
        return;
    }

    int density = maxD_inv + random() % minD_inv;

    for (uint64_t i = 0; i < maxVals; ++i) {
        if (logData && !(i % 100)) {
            printf("\rGenerating matrix. %.0f%% done", 100*((double)i)/maxVals);
            fflush(stdout);
        }
        if (!(i % dest->rows) && (i / dest->rows) < dest->columns) {
            dest->colPtr[(i/dest->rows) + 1] = dest->colPtr[i / dest->rows];
        } 

        if (!(rand() % density)) {
            float nextVal;
            do {
                nextVal = ((float)(rand() % 10000))/100;
            } while(!nextVal);
            vals[dest->valueCount] = nextVal; 
            rowIndices[dest->valueCount++] = i % dest->rows;
            dest->colPtr[(i / dest->rows) + 1]++;
        }
    }

    dest->rowIndices = malloc(dest->valueCount * sizeof(uint64_t));
    dest->values = malloc(dest->valueCount * sizeof(float));

    if (!dest->rowIndices || !dest->values) {
        free(dest->colPtr);
        free(dest->rowIndices);
        free(dest->values);
        errno = ENOMEM;
        return;
    }

    for (uint64_t i = 0; i < dest->valueCount; ++i) {
        if (logData && !(i % 100)) {
            printf("\rStoring data. %.0f%% done.", 100*((double)i)/dest->valueCount);
            fflush(stdout);
        }
        dest->rowIndices[i] = rowIndices[i];
        dest->values[i] = vals[i];
    }
    free(vals);
    free(rowIndices);
}

void generateCSCMatrixTRand(struct cscMatrixTranspose* dest, int minD_inv,
        int maxD_inv) {    

    dest->valueCount = 0;
    uint64_t maxVals;
    bool of = __builtin_umull_overflow(dest->rows, dest->columns, &maxVals);
    if (of) {
        errno = ERANGE;
        return;
    }
    float* vals = malloc(maxVals * sizeof(float));
    uint64_t* rowIndices = malloc(maxVals * sizeof(uint64_t));
    uint64_t* colIndices = malloc(maxVals * sizeof(uint64_t));
    dest->colPtr = calloc(dest->columns+1, sizeof(uint64_t));
    if (!vals || !rowIndices || !dest->colPtr) {
        free(vals);
        free(rowIndices);
        free(colIndices);
        free(dest->colPtr);
        errno = ENOMEM;
        return;
    }

    int density = maxD_inv + random() % minD_inv;

    for (uint64_t i = 0; i < maxVals; ++i) {
        if (logData && !(i % 100)) {
            printf("\rGenerating matrix. %.0f%% done", 100*((double)i)/maxVals);
            fflush(stdout);
        }
        if (!(i % dest->rows) && (i / dest->rows) < dest->columns) {
            dest->colPtr[(i/dest->rows) + 1] = dest->colPtr[i / dest->rows];
        } 

        if (!(rand() % density)) {
            float nextVal;
            do {
                nextVal = ((float)(rand() % 10000))/100;
            } while(!nextVal);
            vals[dest->valueCount] = nextVal; 
            colIndices[dest->valueCount] = i / dest->rows;
            rowIndices[dest->valueCount++] = i % dest->rows;
            dest->colPtr[(i / dest->rows) + 1]++;
        }
    }

    dest->rowIndices = malloc(dest->valueCount * sizeof(uint64_t));
    dest->colIndices = malloc(dest->valueCount * sizeof(uint64_t));
    dest->values = malloc(dest->valueCount * sizeof(float));

    if (!dest->rowIndices || !dest->values || !dest->colIndices) {
        free(dest->colPtr);
        free(dest->rowIndices);
        free(dest->colIndices);
        free(dest->values);
        errno = ENOMEM;
        return;
    }

    for (uint64_t i = 0; i < dest->valueCount; ++i) {
        if (logData && !(i % 100)) {
            printf("\rStoring data. %.0f%% done.", 100*((double)i)/dest->valueCount);
            fflush(stdout);
        }
        dest->rowIndices[i] = rowIndices[i];
        dest->colIndices[i] = colIndices[i];
        dest->values[i] = vals[i];
    }
    free(vals);
    free(rowIndices);
    free(colIndices);
}


/**
 * Stores a row of a CSC matrix into a float*
 *
 * @param m         The desired matrix
 * @param row       The index of the row to store
 * @param result    float* where the row will be stored
 *
 * If row >= m->rows, or the length of result is less than m->columns,
 * this function results in undefined behavior.
 */
static void getFullRow(struct cscMatrix* m, uint64_t row, float* result) {
    for (uint64_t i = 0; i < m->columns; ++i) {
        uint64_t len = m->colPtr[i+1] - m->colPtr[i];
        float col[len];
        uint64_t rowIdxs[len];
        get_col_vector(m->values, m->rowIndices, len, col, rowIdxs, m->colPtr[i]);
        int found = 0;
        for (uint64_t j = 0; j < len; ++j) {
            if (rowIdxs[j] == row) {
                result[i] = col[j];
                found = 1;
                break;
            }       
        }
        if (found) {
            found = 0;
        } else {
            result[i] = 0;
        }
    }
}

void printCSCMatrix(struct cscMatrix* m, int format) {
    char sep = format ? ',' : ' ';
    char rowEnd = format ? '}' : '\n';
    uint64_t n = m->rows;

    if (format) printf ("{");

    for (uint64_t i = 0; i < n; ++i) {
        float row[m->columns];
        getFullRow(m, i, row);
        if (format) printf("{");
        for (uint64_t j = 0; j < m->columns-1; ++j) {
            printf("%.2f%c", row[j], sep);
        }
        printf("%.2f%c", row[m->columns-1], rowEnd);
        if (format && i < n - 1) printf(",");
    }
    if (format) printf("}");
}

void printUint64Vector(uint64_t* arr, uint64_t len) {
    printf("{");
    for (uint64_t i = 0; i < len; ++i) {
        printf("%lu", arr[i]);
        if (i < len-1) {
            printf(", ");
        }
    }
    printf("}\n");
}

void printFloatVector(float* arr, uint64_t len) { 
    printf("{");
    for (uint64_t i = 0; i < len; ++i) {
        printf("%g", arr[i]);
        if (i < len-1) {
            printf(", ");
        }
    }
    printf("}\n");
}

void print_WolframAlpha_random_matrix(int rows, int columns, int density_inv) {
    int nonzero = 0;
    printf("{");
    for(int i = 0; i < rows; ++i) {
        printf("{");
        for(int j = 0; j < columns; ++j) {
            if (rand() % density_inv) {
                printf("0");
            } else {
                printf("%d.%d", (rand() % 100), (rand() % 100));
                nonzero++;
            }
            if (j < columns-1) printf(",");
        }
        printf("}");
        if (i < rows-1) printf(",");
    }
    printf("}\n");
    printf("Density: %.1f%%\n", 100*((float)nonzero)/(rows*columns));
}

void print_random_floats(int n, int prob) {
    int nonzero = 0;
    for (int i = 0; i < n; ++i) {
        if(!(rand() % prob)) {
            printf("%d.%d ", (rand() % 100), (rand() % 100));
            nonzero++;
        } else {
            printf("0 ");
        }
    }
    printf("\n");
    printf("Density: %.1f%%\n", 100*((float) nonzero)/n);
}

void printCSCMatrixMetadata(struct cscMatrix* m) {
    printf("rows: %lu\n", m->rows);
    printf("columns: %lu\n", m->columns);
    printf("valueCount: %lu\n", m->valueCount);
    printf("values: ");
    printFloatVector(m->values, m->valueCount);
    printf("rowIndices: ");
    printUint64Vector(m->rowIndices, m->valueCount);
    printf("colPtr: ");
    printUint64Vector(m->colPtr, m->columns+1);
}
void printCSCMatrixMetadataT(struct cscMatrixTranspose* m) {
    printf("rows: %lu\n", m->rows);
    printf("columns: %lu\n", m->columns);
    printf("valueCount: %lu\n", m->valueCount);
    printf("values: ");
    printFloatVector(m->values, m->valueCount);
    printf("rowIndices: ");
    printUint64Vector(m->rowIndices, m->valueCount);
    printf("colPtr: ");
    printUint64Vector(m->colPtr, m->columns+1);
    printf("colIndices: ");
    printUint64Vector(m->colIndices, m->valueCount);
}

int testVectorPurity(float* vector, uint64_t len) {
    if (!vector) return 1;
    for (uint64_t i = 0; i < len; ++i) {
        if (vector[i] == 0) return 0;
    }
    return 1;
}

int is_valid_csc(struct cscMatrix* m) {
    if (m->rows < 1 || m->columns < 1) return -1;

    for (uint64_t i = 0; i < m->valueCount; ++i) {
        if (m->rowIndices[i] >= m->rows) return -2;
        if (m->values[i] == 0) {
            printf("Zero at index %lu: %f\n", i, m->values[i]);
            return -3;
        }
    }

    for (uint64_t i = 0; i < m->columns; ++i) {
        if (m->colPtr[i+1] < m->colPtr[i]) return -4;
    }
    if (m->colPtr[m->columns] != m->valueCount) return -5;
    return 1;
}

int is_valid_csc_transpose(struct cscMatrixTranspose* m) {
    if (m->rows < 1 || m->columns < 1) return -1;

    for (uint64_t i = 0; i < m->valueCount; ++i) {
        if (m->rowIndices[i] >= m->rows) return -2;
        if (m->colIndices[i] >= m->columns) return -3; 
        if (m->values[i] == 0) {
            printf("Zero at index %lu: %f\n", i, m->values[i]);
            return -4;
        }
        if (i < m->valueCount - 1 && m->colIndices[i+1] < m->colIndices[i]) 
            return -5;
    }

    for (uint64_t i = 0; i < m->columns; ++i) {
        if (m->colPtr[i+1] < m->colPtr[i]) return -6;
    }
    if (m->colPtr[m->columns] != m->valueCount) return -7;
    return 1;
}
