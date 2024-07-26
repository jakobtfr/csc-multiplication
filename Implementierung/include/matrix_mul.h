#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H
#include <stdint.h>

/**
 * Computes A*B with transpose(A) and in-place scalar product.
 *
 * @param a         Matrix A transposed. Passed as struct cscMatrix*  
 * @param b         Matrix B. Passed as struct cscMatrix*
 * @param result    struct cscMatrix* to store result. 
 *                  result->rowIndices, result->colPtr and result->values are 
 *                  stored on the heap.
 */
void matr_mult_csc(const void* a, const void* b, void* result);

/**
 * Computes A*B with transpose(A) and out-of-place scalar product
 *
 * All parameters are identical to those of matr_mult_csc
 */
void matr_mult_csc_V1(const void* a, const void* b, void* result);

/**
 * Computes A*B with no alterations to either matrix (i.e. no transposition)
 * and out-of-place scalar product
 *
 * @param a         Matrix A with no transposition. Passed as struct cscMatrix*  
 * Further parameters are identical to those of matr_mult_csc
 */
void matr_mult_csc_V2(const void* a, const void* b, void* result);

void matr_mult_dense(const void* a, const void* b, void* result, uint64_t a_rows, uint64_t a_cols, uint64_t b_cols);

#endif
