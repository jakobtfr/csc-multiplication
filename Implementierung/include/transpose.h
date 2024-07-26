#ifndef TRANSPOSE_H
#define TRANSPOSE_H
#include <stdint.h>
#include "cs_matrix.h"

int computeColPtr(uint64_t* indices, uint64_t* dest, uint64_t n, uint64_t valueCount);

int transpose(struct cscMatrixTranspose* a, struct cscMatrix* a_t);
/**
 * Copies the contents of a cscMatrixTranspose into a cscMatrix
 *
 * @param a     The cscMatrixTranspose
 * @param a_t   The cscMatrix
 * @return      1 if the operation succeeded, 0 otherwise
 */
int csc_matr_cpy(struct cscMatrixTranspose* a, struct cscMatrix* a_t);

#endif
