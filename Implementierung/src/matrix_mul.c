#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "matrix_mul.h"
#include "cs_matrix.h"

/**
 * Procedure to free all pointer members in a cscMatrix.
 *
 * If any of the pointer members either points to the stack or has already been
 * freed, this function results in undefined behavior.
 *
 * @param result        Matrix whose members should be freed
 */
static void freeResultPtrs(struct cscMatrix* result) {
    free(result->colPtr);
    free(result->rowIndices);
    free(result->values);
}

static void realloc_result(struct cscMatrix* result, uint64_t resultSize){
    if (result->valueCount > 0 && resultSize > result->valueCount) { 
        float* newVals = realloc(result->values, result->valueCount 
                * sizeof(float));
        uint64_t* newInd = realloc(result->rowIndices, result->valueCount
                * sizeof(uint64_t));
        if (!newVals || !newInd) { 
            free(newInd);
            free(newVals);
            freeResultPtrs(result);
            errno = ENOMEM;
            perror("Error reallocating memory for entries of result matrix.");
        }
        result->values = newVals;
        result->rowIndices = newInd;
    } else if (!result->valueCount) {
        free(result->values);
        free(result->rowIndices);
        result->values = 0;
        result->rowIndices = 0;
    }
}

/**
 * Initializes the result matrix of a multiplication with either known values
 * or estimates.
 *
 * The value and rowIndices pointers are initialized to 
 * max(result->rows, result->columns) elements, and the size is returned to allow 
 * dynamic size modifications.
 *
 * @param a             First factor in multiplication 
 * @param b             Second factor in multiplication
 * @param result        Product of first and two factor matrices
 * @param transposed    Boolean integer, nonzero if a is transposed
 * @return              The amount of elements that the pointers point to, or 0 
 *                      if an error ocurred during initialization
 */
static uint64_t initializeResultMatrix(const struct cscMatrix* a, 
        const struct cscMatrix* b, struct cscMatrix* result, int transposed) {

    result->rows = transposed ? a->columns : a->rows;
    result->columns = b->columns;

    uint64_t vals = result->rows > result->columns ? result->rows 
        : result->columns;

    uint64_t maxSize = sizeof(uint64_t) > sizeof(float) ? sizeof(uint64_t) 
        : sizeof(float);
    uint64_t prod;
    uint64_t sum;
    int of_mul = __builtin_umull_overflow(vals, maxSize, &prod);
    int of_add = __builtin_uaddl_overflow(vals, maxSize, &sum);

    if (of_mul || of_add) {
        errno = ERANGE;
        return 0;
    }

    result->valueCount = 0;
    
    result->rowIndices = malloc(vals * sizeof(uint64_t));
    result->values = malloc(vals * sizeof(float));
    result->colPtr = calloc(sum, sizeof(uint64_t));
    if (!result->rowIndices || !result->values || !result->colPtr) {
        freeResultPtrs(result);
        return 0;
    }
    return vals;
}

void matr_mult_csc(const void* a, const void* b, void* result) {
    const struct cscMatrix* csA = a; 
    const struct cscMatrix* csB = b;
    struct cscMatrix* csResult = result;

    errno = 0;
    uint64_t resultSize = initializeResultMatrix(csA, csB, csResult, 1);
    if (errno != 0) {
        perror("Error initializing result matrix members");
        return;
    }
    uint64_t maxSize;
    if (__builtin_umull_overflow(csResult->rows, csResult->columns, &maxSize)) 
        maxSize = UINT64_MAX;

    if (logData) printf("Result matrix members initialized successfully.\n");

    for (uint64_t j = 0; j < csB->columns; ++j) {
        if (logData && csB->columns > 100 && !(j % (csB->columns/100))) {
            printf("\rComputing product of matrices. "
                    "%.0f%% done.", 100*((double) j)/csB->columns);
            fflush(stdout);
        }

        uint64_t bStart = csB->colPtr[j];
        uint64_t bEnd = csB->colPtr[j+1];

        if (bEnd-bStart == 0) {
            if (j < csB->columns-1){
                csResult->colPtr[j+2] = csResult->colPtr[j+1];
            }
            continue;
        }

        for (uint64_t i = 0; i < csA->columns; ++i) {
            // Compute the scalar product between the next A column and the current
            // B column in-place
            float entry = scalar_prod_in_place(csA->values, csB->values, 
                    csA->rowIndices, csB->rowIndices, csA->colPtr[i], 
                    csA->colPtr[i+1], bStart, bEnd);

            if (cmp_float_eq(entry, 0)) continue; 

            // Increase the memory for values and row indices if necessary
            if (csResult->valueCount >= resultSize) {
                resultSize = extend_vector(&csResult->values, 
                        &csResult->rowIndices, csResult->valueCount, maxSize);
                if (!resultSize) {
                    perror("Error storing result values.");
                    freeResultPtrs(csResult);
                    return;
                }
            }

            csResult->values[csResult->valueCount] = entry; 
            csResult->rowIndices[csResult->valueCount++] = i;
            csResult->colPtr[j+1]++;
        }       
        // Set value of next-next column pointer to value of next one
        if (j < csB->columns-1){
            csResult->colPtr[j+2] = csResult->colPtr[j+1];
        }
    }

    if (logData) printf("\rProduct of matrices computed successfully.\n");

    errno = 0;
    // realloc result values and rowIndices to valueCount
    realloc_result(csResult, resultSize);
}

void matr_mult_csc_V1(const void* a, const void* b, void* result) {
    const struct cscMatrix* csA = a; 
    const struct cscMatrix* csB = b;
    struct cscMatrix* csResult = result;

    errno = 0;
    uint64_t resultSize = initializeResultMatrix(csA, csB, csResult, 1);
    if (errno != 0) {
        perror("Error initializing result matrix members");
        return;
    }
    uint64_t maxSize;
    if (__builtin_umull_overflow(csResult->rows, csResult->columns, &maxSize)) 
        maxSize = UINT64_MAX;

    if (logData) printf("Result matrix members initialized successfully.\n");

    for (uint64_t j = 0; j < csB->columns; ++j) {
        if (logData && csB->columns > 100 && !(j % (csB->columns/100))) {
            printf("\rComputing product of matrices. "
                    "%.0f%% done.", 100*((double) j)/csB->columns);
            fflush(stdout);
        }
        // get the next B column
        uint64_t start = csB->colPtr[j];
        uint64_t end = csB->colPtr[j+1];
        uint64_t bLen = end-start;
        if (bLen == 0) {
            if (j < csB->columns-1){
                csResult->colPtr[j+2] = csResult->colPtr[j+1];
            }
            continue;
        }

        float* bVec = malloc(bLen * sizeof(float));
        uint64_t* bInd = malloc(bLen * sizeof(uint64_t));
        if (!bVec || !bInd) {
            errno = ENOMEM;
            free(bVec);
            free(bInd);
            freeResultPtrs(csResult);
            fprintf(stderr, 
                    "\rError allocating memory for column %lu of matrix B.", j);
            return;
        }

        get_col_vector(csB->values, csB->rowIndices, bLen, bVec, bInd, start);

        for (uint64_t i = 0; i < csA->columns; ++i) {
            start = csA->colPtr[i];
            end = csA->colPtr[i+1];
            uint64_t aLen = end-start;
            if (!aLen) continue;

            // get the next A column
            float* aVec = malloc(aLen * sizeof(float));
            uint64_t* aInd = malloc(aLen * sizeof(uint64_t));
            if (!aVec || !aInd) {
                errno = ENOMEM;
                free(aVec);
                free(aInd);
                free(bVec);
                free(bInd);
                freeResultPtrs(csResult);
                fprintf(stderr, 
                    "\rError allocating memory for row %lu of matrix A.\n", i);
                return;
            }

            get_col_vector(csA->values, csA->rowIndices, aLen, aVec, aInd, 
                    start);

            float entry = scalar_prod(aVec, bVec, aInd, bInd, aLen, bLen);

            free(aVec);
            free(aInd);

            if (cmp_float_eq(entry, 0)) continue; 

            // Increase the memory for values and row indices if necessary
            if (csResult->valueCount >= resultSize) {
                resultSize = extend_vector(&csResult->values, &csResult->rowIndices,
                                           csResult->valueCount, maxSize);
                if (!resultSize) {
                    perror("Error storing result values.");
                    free(bVec);
                    free(bInd);
                    freeResultPtrs(csResult);
                    return;
                }
            }

            csResult->values[csResult->valueCount] = entry; 
            csResult->rowIndices[csResult->valueCount++] = i;
            csResult->colPtr[j+1]++;
        }       
        
        // Set value of next-next column pointer to value of next one
        if (j < csB->columns-1){
            csResult->colPtr[j+2] = csResult->colPtr[j+1];
        }
        free(bVec);
        free(bInd);
    }
    if (logData) printf("\rProduct of matrices computed successfully.\n");

    errno = 0;
    // realloc result values and rowIndices to valueCount
    realloc_result(csResult, resultSize);
    
}


void matr_mult_csc_V2(const void* a, const void* b, void* result) {
    const struct cscMatrix* csA = a; 
    const struct cscMatrix* csB = b;
    struct cscMatrix* csResult = result;

    errno = 0;
    uint64_t resultSize = initializeResultMatrix(csA, csB, csResult, 0);
    if (errno != 0) {
        perror("Error initializing result matrix members");
        return;
    }
    uint64_t maxSize;
    if (__builtin_umull_overflow(csResult->rows, csResult->columns, &maxSize)) 
        maxSize = UINT64_MAX;

    if (logData) printf("Result matrix members initialized successfully.\n");

    for (uint64_t j = 0; j < csB->columns; ++j) {
        // show progress if l option is set. Print only 1% of iterations to 
        // reduce overhead
        if (logData && !(j % 100)) {
            printf("\rComputing product of matrices. "
                    "%.0f%% done.", 100*((double) j)/csB->columns);
            fflush(stdout);
        }                                                                        
        // get the next B column
        uint64_t start = csB->colPtr[j];
        uint64_t end = csB->colPtr[j+1];
        uint64_t bLen = end-start;
        if (bLen == 0) {
            if (j < csB->columns-1){
                csResult->colPtr[j+2] = csResult->colPtr[j+1];
            }
            continue;
        }

        float* bVec = malloc(bLen * sizeof(float));
        uint64_t* bInd = malloc(bLen * sizeof(uint64_t));
        if (!bVec || !bInd) {
            errno = ENOMEM;
            free(bVec);
            free(bInd);
            freeResultPtrs(csResult);
            fprintf(stderr, 
                    "\rError allocating memory for column %lu of matrix B.", j);
            return;
        }
                                                                                               
        get_col_vector(csB->values, csB->rowIndices, bLen, bVec, bInd, start);
                                                                               
        for (uint64_t i = 0; i < csA->rows; ++i) {

            errno = 0;
            struct cscRow row = get_row_vector(csA, i);
            if (errno) {
                perror("Error getting row vector from a");
                free(bVec);
                free(bInd);
                freeResultPtrs(csResult);
                return;
            }
            // printFloatVector(row.values, row.valueCount);
            // printUint64Vector(row.colIndices, row.valueCount);

            float entry = scalar_prod(row.values, bVec, row.colIndices, bInd,
                                      row.valueCount, bLen);
            free(row.colIndices);
            free(row.values);

            if (cmp_float_eq(entry, 0)) continue; 

            // Increase the memory for values and row indices if necessary
            if (csResult->valueCount >= resultSize) {
                resultSize = extend_vector(&csResult->values, &csResult->rowIndices,
                                           csResult->valueCount, maxSize);
                if (!resultSize) {
                    perror("Error storing result values.");
                    free(bVec);
                    free(bInd);
                    freeResultPtrs(csResult);
                    return;
                }
            }

            csResult->values[csResult->valueCount] = entry; 
            csResult->rowIndices[csResult->valueCount++] = i;
            csResult->colPtr[j+1]++;
        }       
        
        // Set value of next-next column pointer to value of next one
        if (j < csB->columns-1){
            csResult->colPtr[j+2] = csResult->colPtr[j+1];
        }
        free(bVec);
        free(bInd);
    }
    if (logData) printf("\rProduct of matrices computed successfully.\n");

    errno = 0;
    // realloc result values and rowIndices to valueCount
    realloc_result(csResult, resultSize);
}

