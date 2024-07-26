#ifndef CSC_IO_H
#define CSC_IO_H

#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include "cs_matrix.h"

extern const char* usage_msg;

extern const char* help_msg;

extern const char* shortopts;

extern const struct option longopts[];

/**
 * Prints content of usage_msg.
 *
 *@param progname          Name of the program. 
 */
void print_usage(const char* progname);

/**
 * Prints content of usage_msg and help_msg.
 *
 * @param progname          Name of the program. 
 */
void print_help(const char* progname);

/**
 * Frees the memory allocated for the matrix and its attributes.
 *
 * @param matrix          CSCMatrix to be freed
 */
void free_csc_matrix(struct cscMatrix* matrix);

/**
 * Frees the memory allocated for the matrix and its pointer members.
 *
 * @param matrix          CSCTransposeMatrix to be freed
 */
void free_transpose_matrix(struct cscMatrixTranspose* matrix);

/**
 * Converts a char into int
 * 
 * @param c             Char value of the number to convert. Passed as char*
 * @param i             Output parameter. Passed as int*
 */
int convert_int(char* c, int* i);

/**
 * Converts a char into unsigned int
 * 
 * @param c             Char value of the number to convert. Passed as char*
 * @param ui            Output parameter. Passed as unsigned int*
 */
int convert_unsigned(char* c, unsigned int* ui);

/**
 * Calculates the column indices from a csc Matrix.
 *
 * @param matrix         Matrix from which the indexes are to be calculated.
 */
void get_col_indices(struct cscMatrixTranspose* matrix, uint64_t* coolIndex);


/**
 * Parses file into cscMatrix struct
 *
 * @param filename_a         Filename of a Matrix. 
 * @param filename_b         Filename of a Matrix. 
 * @param matrixA            Matrix in which the values of the documents to be saved.
 * @param matrixB            Matrix in which the values of the documents to be saved.
 * @return                   Returns if one of the matrices is only made up of zeros.
 */
int parse_csc_file(const char* filename_a, const char* filename_b, 
        struct cscMatrixTranspose* matrixA, struct cscMatrix* matrixB);

int parse_csc_file_V2(const char* filename_a, const char* filename_b, 
        struct cscMatrix* matrixA, struct cscMatrix* matrixB);

/**
 * Parses result into the output file.
 *
 * @param result_matrix       Result of the multiplication. Passed as const cscMatrix*.  
 * @param output_file         Filename of the output file. Passed as const char*.  
 */
void result_to_file(struct cscMatrix* result, const char* output_file);

/**
 * Creates a random matrix and writes it in a document in the format csc.
 *
 * @param maxSize             Sets the maximum size for rows and columns for the matrix.  
 * @param minSize             Sets the minimum size for rows and columns for the matrix.
 */
void create_random_matrix_doc(uint64_t minSize, uint64_t maxSize, const char* filenameA, const char* filenameB);

/**
 * Writes a matrix without values other than zero in csc format in a file.
 *
 * @param rows          Rows of the matrix.
 * @param columns       Columns of the matrix.
 */
void print_empty_matrix(uint64_t rows, uint64_t columns, const char* output_file);

#endif
