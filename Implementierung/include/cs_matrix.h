#ifndef CS_MATRIX_H
#define CS_MATRIX_H

#include <stdint.h>

extern int logData;

/**
 * @class cscMatrix
 *
 * Compressed Sparse Column Matrix
 *
 * @member rows         Amount of rows
 * @member columns      Amount of columns
 * @member valueCount   Amount of entries in values and rowIndices
 * @member rowIndices   Row index of each value
 * @member colPtr       Array that splits values based on which column they are in.
 *                      It has m+1 elements in a n*m matrix by definition.
 * @member values       Nonzero values in the matrix
 */
struct cscMatrix {
    uint64_t rows;
    uint64_t columns;
    uint64_t valueCount;
    float* values;
    uint64_t* rowIndices;
    uint64_t* colPtr;
};

/**
 * @class cscMatrixTranspose
 *
 * CSC Matrix for transposing
 *
 * @member rows         Amount of rows
 * @member columns      Amount of columns
 * @member valueCount   Amount of entries in values and rowIndices
 * @member values       Nonzero values in the matrix
 * @member rowIndices   Row index of each value
 * @member colIndices   Column index of each  value
 * @member colPtr       Array that splits values based on which column they are in.
 *                      It has m+1 elements in a n*m matrix by definition.
 */
struct cscMatrixTranspose {
    uint64_t rows;
    uint64_t columns;
    uint64_t valueCount;
    float* values;
    uint64_t* rowIndices;
    uint64_t* colIndices;
    uint64_t* colPtr;
};

/**
 * @class cscRow
 *
 * Row vector in a CSC Matrix
 *
 * @member valueCount   Amount of entries in values and colIndices
 * @member colIndices   Column indices of each value, analogous to a column
 *                      vector's row indices
 * @member values       Nonzero values in the row 
 *
 */
struct cscRow {
    uint64_t valueCount;
    uint64_t* colIndices;
    float* values;
};

/**
 * Gets the column vector from values specified by indexes start and start+1
 * in the vector's colPtr, as well as the row indices of the
 * vector's elements.
 *
 * @param values        The CSC Matrix values
 * @param indices       Row indices of values
 * @param len           Length of vector and vecIndices
 * @param vector        Array to write the values to 
 * @param vecindices    Array to write the row indices to
 * @param start         Starting index in values
 */
void get_col_vector(float* values, uint64_t* indices, uint64_t len,
                    float* vector, uint64_t* subIndices, uint64_t start);
/**
 * Increases the size of the memory pointed to by a values and indices pointer
 * pair
 * The new size will be min(2*n, maxSize)
 *
 * @param values        The values pointer to extend
 * @param indices       The indices pointer to extend
 * @param current       Current amount of bytes allocated to each pointer
 * @param max           Maximum amount of bytes computable to allocate.
 *                      Equivalent to min(m->rows * m->columns, UINT64_MAX)
 * @return              The new allocated element count if succesful, 0 otherwise
 */
uint64_t extend_vector(float** values, uint64_t** indices, uint64_t current,
                       uint64_t max);

/**
 * Gets a row vector from a struct cscMatrix and its column indices (equivalent
 * to a column vector from the transposed matrix). 
 *
 * If an error occurs, errno is set and the returned struct's pointer members
 * are null
 *
 * @param m             The CSC Matrix
 * @param index         The index of the row to get
 * @param vector        Array to write the values to
 * @param vecIndices    Array to write the column indices to
 * @return              The extracted row. The struct's pointer members are 
 *                      stored on the heap
 */
struct cscRow get_row_vector(const struct cscMatrix* m, uint64_t index);

/**
 * Computes the scalar product between two vectors, skipping indices where either
 * vector has a zero.
 *
 * @param aVec      The first vector
 * @param bVec      The second vector
 * @param aInd      Index array of aVec's nonzero values
 * @param bInd      Index array of bVec's nonzero values
 * @param aLen      Length of a
 * @param bLen      Length of b
 * @return          The scalar product
 */
float scalar_prod(const float* aVec, const float* bVec, const uint64_t* aInd,
                  const uint64_t* bInd, uint64_t aLen, uint64_t bLen);

/**
 * Computes the scalar product between two columns of two CSC matrices.
 * The columns are given through column pointer slicing.
 *
 * @param a         The first matrix
 * @param b         The second matrix
 * @param aStart    The starting index of the column of a 
 * @param aEnd      The ending index of the column of a 
 * @param bStart    The starting index of the column of b 
 * @param bEnd      The ending index of the column of b 
 * @return          The scalar product between the columns
 *
 * The starting and ending indices of the ith column of a CSC matrix are given
 * by colPtr[i] and colPtr[i+1] (exclusive).
 *
 */
float scalar_prod_in_place(const float* aVec, const float* bVec, 
        const uint64_t* aInd, const uint64_t* bInd, uint64_t aStart, 
        uint64_t aEnd, uint64_t bStart, uint64_t bEnd);
/**
 * Compares two struct cscMatrix a and b for equality. a = b iff 
 * every member of a is equal to the equivalent member in b.
 *
 * @param a     First matrix
 * @param b     Second matrix
 * @return      1 if a = b, 0 otherwise
 */
int cmp_csc_eq(struct cscMatrix* a, struct cscMatrix* b);

/**
 * Compares two uint64_t vectors for equality.
 *
 * @param a     The first vector
 * @param b     The second vector 
 * @param n     The length of a and b
 * @return      1 if a[i] == b[i] for every i s.t. 0 <= i < n, 0 otherwise
 *
 * If n > a_len or n > b_len this function results in undefined behavior.
 *
 * If n < a_len or n < b_len the function only compares the vectors from 
 * index 0 to n-1, without checking the intervals from n to a_len and/or 
 * from n to b_len.
 *
 */
int cmp_uint64_t_vec_eq(uint64_t* a, uint64_t* b, uint64_t n);

/**
 * Compares two float vectors for equality with a margin of error as defined 
 * in the function cmp_float_eq.
 *
 * @param a     The first vector
 * @param b     The second vector
 * @param n     The length of a and b
 * @return      1 if cmp_float_eq(a[i], b[i]) returns 1 for every i s.t. 0 <= i < n,
 *              0 otherwise.
 *              See the documentation of cmp_uint64_t_vec_eq for further information.
 *         
 */
int cmp_float_vec_eq(float* a, float* b, uint64_t n);

/**
 * Generates a n*n CSC identity matrix
 *
 * @param id    Struct to store matrix
 * @param n     Row and column count of matrix
 */
int generate_csc_id(struct cscMatrix* id, uint64_t n);

/**
 * Generates a cscMatrix with the given parameters. Refer to the documentation
 * of struct cscMatrix for information about the struct members.
 */
void generate_csc_matr(struct cscMatrix* dest, uint64_t rows, uint64_t columns,
                       uint64_t valueCount, uint64_t* rowIndices, uint64_t* colPtr,
                       float* values);


/**
 * Generates a cscMatrixTranspose with the given parameters. Refer to the documentation
 * of struct cscMatrixTranspose for information about the struct members.
 */
void generateCSCMatrixTranspose(struct cscMatrixTranspose* dest, uint64_t rows, uint64_t columns,
                       uint64_t valueCount, uint64_t* rowIndices, uint64_t* colIndices, uint64_t* colPtr,
                       float* values);

/**
 * Calculates the column indices for a given array of column pointers.
 * @param colPtr    The given column pointers
 * @param numCols   Number of columns
 * @param colInd    The array where the computed column indices are stored
 */
void calculateColumnIndices(uint64_t *colPtr, uint64_t numCols, uint64_t *colInd);


/**
 * Prints the given matrix as a string, either with every row separated by a 
 * newline and every column by a whitespace character, or in Wolfram language
 * plain text.
 *
 * @param m         Matrix to print.
 * @param format    Format to print the matrix in. If set to 0, whitespaces and 
 *                  newlines are used, otherwise, it will be printed in Wolfram
 *                  language format (does not end in newline). 
 *
 *                  Refer to the documentation of 
 *                  print_WolframAlpha_random_matrix for further information on
 *                  the Wolfram language format.
 */
void printCSCMatrix(struct cscMatrix* m, int format);

/**
 * Prints the contents of a uint64[]
 *
 * @param arr   The array to print
 * @param len   The length of the array
 *
 * If len is greater than the length of arr, this function results in undefined
 * behavior.
 * If len is less than the length of arr, only the contents from index 0 to 
 * len-1 in arr will be printed.
 */
void printUint64Vector(uint64_t* arr, uint64_t len);


void printFloatVector(float* arr, uint64_t len);

/**
 * Compares floats for equality with a fixed margin of error.
 * Algorithm from https://floating-point-gui.de/errors/comparison/
 *
 * @param a     The first float to compare
 * @param b     The second float to compare
 * @return      1 if a = b with margin of error of approximately 0.001%, 
 *              0 otherwise
 */
int cmp_float_eq(float a, float b);

/**
 * Prints a randomized matrix in plain Wolfram language text to the console. 
 * The matrix's density is printed underneath. 
 *
 * The printed matrix can be pasted directly into WolframAlpha with no further 
 * modifications.
 *
 * @param rows          The amount of rows in the matrix
 * @param columns       The amount of columns in the matrix
 * @param density_inv   The reciprocal of the desired density of the matrix.
 *                      Since the standard library's random number generation 
 *                      does not adequately simulate randomness, multiple function
 *                      calls may be necessary to reach the desired density.
 */
void print_WolframAlpha_random_matrix(int rows, int columns, int density_inv);


/**
 * Prints a string of random floats, most of them being zeros.
 *
 * @param n         The amount of floats to print
 * @param prob      The reciprocal of the percentage of nonzero numbers 
 *                  See density_inv in print_WolframAlpha_random_matrix for more
 *                  information.
 */
void print_random_floats(int n, int prob);

/**
 * Generates a random CSC Matrix with expected density between given ranges. 
 * Because of the standard library's poor random number generation, the density
 * can be below or above these boundaries.
 *
 * The matrix's dimensions must be passed within the struct; i.e., m->rows and
 * m->columns must be initialized. It is recommended to also use some form of
 * random number generation for this.
 *
 * All pointer struct members are stored on the heap.
 * If memory allocation causes an error, errno will be set to ENOMEM.
 *
 * @param dest      Pointer to the struct where the matrix will be stored.
 * @param minD_inv  The reciprocal of the minimum desired density
 * @param maxD_inv  The reciprocal of the maximum desired density
 */
void generate_csc_matr_rand(struct cscMatrix* dest, int minD_inv, int maxD_inv);

void generateCSCMatrixTRand(struct cscMatrixTranspose* dest, int minD_inv,
        int maxD_inv);    

/**
 * Prints the members of a struct cscMatrix
 *
 * @param m     Matrix to print
 */
void printCSCMatrixMetadata(struct cscMatrix* m);
void printCSCMatrixMetadataT(struct cscMatrixTranspose* m);

/**
 * Tests the purity of a vector. A vector is considered pure if it contains no
 * zeros.
 *
 * @param vector    The vector to test
 * @param len       The vector's length
 * @return          1 if the vector is pure, 0 otherwise
 */
int testVectorPurity(float* vector, uint64_t len);

/**
 * Checks if the metadata in a struct cscMatrix corresponds to a valid CSC matrix
 *
 * Since array lengths are unknown by the compiler, the function tests for valid
 * array lengths by performing reads from 0 to the array length the array is 
 * supposed to have (e.g. valueCount for rowIndices or columns+1 for colPtr).
 *
 * To allow proper testing, valgrind or a similar tool should be used, since 
 * it will show when an invalid read (i.e. array too short) is performed
 * (this can result in a false positive if valgrind isn't used). If the arrays
 * happen to be larger than the value count, the function can return a false
 * positive.
 *
 * @param m     The matrix to test
 * @return      1 if the matrix is valid, 
 *              otherwise:
 *              -1: invalid dimensions
 *              -2: row indices out of range
 *              -3: zero in value ptr
 *              -4: colPtr not sorted
 *              -5: last index of colPtr is not valueCount
 */
int is_valid_csc(struct cscMatrix* m);

/**
 * Checks if the metadata in a struct cscMatrixTranspose corresponds to a valid 
 * CSC matrix
 *
 * See the documentation of is_valid_csc for further information about the behavior
 * of this function.
 *
 * @param m     The matrix to test
 * @return      1 if the matrix is valid, 
 *              otherwise,
 *              -1: invalid dimensions
 *              -2: row indices out of range
 *              -3: column indices out of range
 *              -4: zero in value ptr
 *              -5: column indices not sorted
 *              -6: colPtr not sorted
 *              -7: last index of colPtr is not valueCount
 */
int is_valid_csc_transpose(struct cscMatrixTranspose* m);

#endif
