#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>

#include "csc_io.h"
#include "cs_matrix.h"

const char* usage_msg =
    "Usage: %s [commands] <arg>\n"
    "   or: %s -h, --help             Show help message and exits.\n"
    "The use of commands is optional. If you choose not to use a command a"
    " standard value will be used for the execution of the program.\n";


const char* help_msg =
    "Help Message\n"
    "Commands with mandatory arguments:\n"
    "  -V <Number>          Specify the version of the function.\n"
    "  -a <Filename>        Specify file containing Matrix a.\n"
    "  -b <Filename>        Specify file containing Matrix b.\n"
    "  -o <Filename>        Specify the output file .\n"
    "Commands with optional arguments:\n"
    "  -B<N>                Benchmarking mode. Logs the execution time of the "
                            "program to the console, as well as the duration of"
                            " different operations.\n" 
    "                       N specifies the amount of times to perform the multiplication.\n"
    "Commands without arguments:\n"
    "  -h, --help           Display this help message and exits.\n"
    "  -l                   Prints messages to the console indicating the "
                            "progress of the program.\n"
    "                       This option adds significant overhead and "
                            "increases its runtime by a nonnegligible amount.\n"
    "                       Not recommended for benchmarking mode.\n"
    "  -r                   Generates a random pair of matrices A and B with "
                            "dimensions n*m and m*k respectively,\n"
    "                       for n,m,k in between 512 and 1024. "
                            "The resulting matrices are written to the files\n"
    "                       randomMatrixA.txt and randomMatrixB.txt.\n";

const char* shortopts = "V:a:b:o:B::hlr";

const struct option longopts[] = {
    {"help", no_argument, 0, 'h'},
    {0,0,0,0}
};


void print_usage(const char* progname) {
    fprintf(stderr, usage_msg, progname, progname);
}


void create_random_matrix_doc(uint64_t minSize, uint64_t maxSize, const char* filenameA, const char* filenameB) {
    struct cscMatrix matrixA = {0};
    struct cscMatrix matrixB = {0};

    uint64_t diff = maxSize - minSize + 1;

    uint64_t aRows = minSize + rand() % diff;
    uint64_t aColbRows = minSize + rand() % diff;
    uint64_t bCols = minSize + rand() % diff;

    matrixA.columns = matrixB.rows = aColbRows;
    matrixA.rows = aRows;
    matrixB.columns = bCols;

    if (logData) {
        printf("Generating matrix A...");
        fflush(stdout);
    }
    errno = 0;
    generate_csc_matr_rand(&matrixA, 10, 3);
    if (errno) {
        perror("\rError while generating matrix A");
        return;
    }
    if (logData) {
        printf("\rMatrix A generated successfully.\n");
        printf("Generating matrix B...");
        fflush(stdout);
    }
    errno = 0;
    generate_csc_matr_rand(&matrixB, 10, 3);
    if (errno) {
        perror("\rError while generating matrix B");
        return;
    }
    if (logData) {
        printf("\rMatrix B generated successfully.\n");
    }

    if (logData) {
        printf ("Writing matrix A into input file...");
        fflush(stdout);
    }
    errno = 0;
    result_to_file(&matrixA, filenameA);
    if (errno) return;
    if (logData) {
        printf("\rMatrix A written successfully.     \n");
        printf ("Writing matrix B into input file...");
        fflush(stdout);
    }
    result_to_file(&matrixB, filenameB);
    if (errno) return;
    if (logData) printf("\rMatrix B written successfully.     \n");
    free(matrixA.colPtr);
    free(matrixA.values);
    free(matrixA.rowIndices);
    free(matrixB.colPtr);
    free(matrixB.values);
    free(matrixB.rowIndices);
}


void print_help(const char* progname) {
    print_usage(progname);
    fprintf(stderr, "\n%s", help_msg);
}


void free_csc_matrix(struct cscMatrix* matrix){
    free(matrix->colPtr);
    free(matrix->values);
    free(matrix->rowIndices);
    free(matrix);
}


void free_transpose_matrix(struct cscMatrixTranspose* matrix){
    free(matrix->values);
    free(matrix->rowIndices);
    free(matrix->colPtr);
    free(matrix->colIndices);
    free(matrix);
}


int convert_int(char* c, int* i) {
    errno = 0;
    char* endptr;

    long val = strtol(c, &endptr, 10);

    if (endptr == c || *endptr != '\0') {
        fprintf(stderr, "Invalid number: %s could not be converted to int\n", c);
        return 1;
    }

    else if (errno == ERANGE || val < INT_MIN || val > INT_MAX) {
        fprintf(stderr, "Invalid number: %s over- or underflows int\n", c);
        return 1;
    }

    *i = (int)val;
    return 0;
}


int convert_unsigned(char* c, unsigned int* ui) {
    errno = 0;
    char* endptr;

    unsigned long val = strtoul(c, &endptr, 10);

    if (endptr == c || *endptr != '\0') {
        fprintf(stderr, "Invalid number: %s could not be converted to unsigned int\n", c);
        return 1;
    }

    else if (errno == ERANGE || val > UINT_MAX) {
        fprintf(stderr, "Invalid number: %s overflows unsigned int\n", c);
        return 1;
    }

    *ui = (unsigned int)val;
    return 0;
}


/**
 * Counts the values in a line of a document
 *
 * @param line      Line to parse
 */
static uint64_t count_values(char* line) {
    if (line[0] == '\n' || line[0] == '\0') {
        return 0;
    }
    uint64_t count = 1;  // There's at least one token if the line is not empty
    for (char* c = line; *c != '\0'; c++) {
        if (*c == ',') {
            count++;
        }
    }
    return count;
}


/**
 * Parses the values of a line and return a float array pointer of the parsed values.
 *
 * @param array       Array where the values are to be stored. Passed as float*.  
 * @param num         Amount of elements in the line. Passed as uint64_t.  
 * @param line        String line to be parsed. Passed as char*.  
 */
static void line_parsing_float(float* array, uint64_t num, char* line) {
    char* value = strtok(line, ",");
    char *endptr;
    float valAsFloat;
    for (uint64_t i = 0; i < num; i++) {
        valAsFloat = strtof(value, &endptr);
        if (valAsFloat == 0)
        {   
            errno = EINVAL;
            perror("Wrong format. The values can not contain a zero value");
            return;
        }
        array[i] = valAsFloat;
        value = strtok(NULL, ",");
    }
}

/**
 * Parses the values of a line and return a uint64_t array pointer of the parsed values.
 *
 * @param array       Array where the values are to be stored. Passed as uint64_t*.  
 * @param num         Amount of elements in the line. Passed as uint64_t.  
 * @param line        String line to be parsed. Passed as char*.  
 * @return            Returns if all colPtr are zero.
 */
static int line_parsing_uint(uint64_t* array, uint64_t num, char* line, char* type, uint64_t rows) {
    char* value = strtok(line, ",");
    char *endptr;
    uint64_t valAsUnit;
    uint64_t allZeros = 1;
    for (uint64_t i = 0; i < num; i++) {
        valAsUnit = strtoull(value, &endptr, 10);
        if (!strcmp(type, "row") && valAsUnit > rows - 1){
            errno = EINVAL;
            fprintf(stderr, "Wrong format. The row indices can not contain a value bigger than: %"PRIu64"", rows - 1);
            return 0;
        }else if (!strcmp(type, "column")) {
            if (valAsUnit != 0){
                allZeros = 0;
            }
        }
        array[i] = valAsUnit;
        value = strtok(NULL, ",");
    }
    return allZeros;
}

/**
 * Calculates the column indices from a csc Matrix.
 *
 * @param matrix         Matrix from which the indexes are to be calculated.
 * @param colIndices     Column indices of the matrix.
 */
void get_col_indices(struct cscMatrixTranspose* matrix, uint64_t* colIndices) {
    for (uint64_t i = 0; i < matrix->columns; ++i) {
        uint64_t start = matrix->colPtr[i];
        uint64_t len = matrix->colPtr[i+1] - matrix->colPtr[i];
        for (uint64_t j = start; j < start+len; ++j) {
            colIndices[j] = i;
        }
    }
}

/**
 * Parses the dimensions into the matrices
 *
 * @param aRows         Rows of MatrixA.
 * @param aCols         Columns of MatrixB.
 * @param bRows         Rows of MatrixB.
 * @param bCols         Columns of MatrixB.
 * @param fileA         File of a Matrix.
 * @param fileB         File of a Matrix.
 */
static void read_dimensions(uint64_t* aRows, uint64_t* aCols, uint64_t* bRows, uint64_t*bCols, FILE* fileA, FILE* fileB) {
    char lineA[1024];
    char lineB[1024];

    if (fgets(lineA, sizeof(lineA), fileA) && fgets(lineB, sizeof(lineB), fileB)) {
        sscanf(lineA, "%lu,%lu", aRows, aCols);
        sscanf(lineB, "%lu,%lu", bRows, bCols);
        
        // Edge cases
        if (*aCols != *bRows) {
            fprintf(stderr, "Dimension mismatch: cannot multiply %lu by %lu "
                    "matrix by a %lu by %lu matrix.\n",*aRows, *aCols, *bRows, *bCols);
            errno = EINVAL;
			return;
        } else if (*aCols <= 0 || *aRows <= 0 || *bCols <= 0 || *bRows <= 0) {
            fprintf(stderr, "Impossible dimension: cannot multiply %lu by %lu "
                    "matrix by a %lu by %lu matrix.\n", *aRows, *aCols, *bRows, *bCols);
            errno = EINVAL;
			return;
        }
    }

}

/**
 * Parses the values into the matrices
 *
 * @param aValC         Value count of MatrixA.
 * @param bValC         Value count of MatrixB.
 * @param aVals         Values of MatrixA.
 * @param bVals         Values of MatrixB.
 * @param fileA         File of a Matrix. 
 * @param fileB         File of a Matrix.
 */
static void read_values(uint64_t* aValC, uint64_t* bValC, float** aVals,
        float** bVals, FILE* fileA, FILE* fileB) {
    char* lineA = NULL;
    char* lineB = NULL;
    size_t lenA = 0;
    size_t lenB = 0;

    getline(&lineA, &lenA, fileA);
    getline(&lineB, &lenB, fileB);

    if (lineA && lineB) { 
        // Count amount of values in each matrix
        *aValC = count_values(lineA);
        *bValC = count_values(lineB);
        if (*aValC == 0 || *bValC == 0)
        {
            *aVals = 0;
            *bVals = 0;
        }else {
            // Allocate memory for the values
            *aVals = malloc(*aValC * sizeof(float));
            *bVals = malloc(*bValC * sizeof(float));
            if (!*aVals || !*bVals) {
                free(*aVals);
                free(*bVals);
                free(lineA);
                free(lineB);
                errno = ENOMEM;
                return;
            }

            // Store values
            line_parsing_float(*aVals, *aValC, lineA);
            if(errno) return;
            line_parsing_float(*bVals, *bValC, lineB);
            if(errno) return;
        }
        free(lineA);
        free(lineB);
    } else {
        perror("An error occurred during parsing");
        free(lineA);
        free(lineB);
    }

}

/**
 * Parses the row indices into the matrices
 *
 * @param aValC         Value count of MatrixA. 
 * @param aRows         Rows of MatrixA.
 * @param aVals         Values of MatrixA.
 * @param aRowIds       Row indices of MatrixA.
 * @param bValC         Value count of MatrixB.
 * @param bRows         Rows of MatrixB.
 * @param bVals         Values of MatrixB.
 * @param bRowIds       Row indices of MatrixB.
 * @param fileA         File of a Matrix.
 * @param fileB         File of a Matrix.
 */
static void read_row_indices(uint64_t aValC, uint64_t aRows, float** aVals, 
        uint64_t** aRowIds, uint64_t bValC, uint64_t bRows, float** bVals, 
        uint64_t** bRowIds, FILE* fileA, FILE* fileB) {
    char* lineA = NULL;
    char* lineB = NULL;
    size_t lenA = 0;
    size_t lenB = 0;

    getline(&lineA, &lenA, fileA);
    getline(&lineB, &lenB, fileB);

    if (lineA && lineB) {
        if (count_values(lineA) != aValC || count_values(lineB) != bValC)
        {
            errno = EINVAL;
            perror("Wrong format. The number of indices and values do no match");
            free(*aVals);
            free(*bVals);
            free(*aRowIds);
            free(bRowIds);
            return;
        }
        if (aValC == 0 || bValC == 0) {
            *aRowIds = 0;
            bRowIds = 0;
        } else {
            // Allocate memory for the row indices
            *aRowIds = malloc(aValC * sizeof(uint64_t));
            *bRowIds = malloc(bValC * sizeof(uint64_t));
            if (!*aRowIds || !*bRowIds) {
                free(*aVals);
                free(*bVals);
                free(*aRowIds);
                free(*bRowIds);
                free(lineA);
                free(lineB);
                errno = ENOMEM;
                return;
            }

            // Store row indices
            line_parsing_uint(*aRowIds, aValC, lineA, "row", aRows);
            line_parsing_uint(*bRowIds, bValC, lineB, "row", bRows);
        }
        free(lineA);
        free(lineB);
    } else {
        perror("An error occurred during parsing");
        free(lineA);
        free(lineB);
    }
}

/**
 * Parses the column pointers into the matrices
 *
 * @param aCols         Columns of MatrixA. 
 * @param aColPtr       Column Pointers of MatrixA.
 * @param bCols         Columns of MatrixB. 
 * @param bColPtr       Column Pointers of MatrixB.
 * @param fileA         File of a Matrix.
 * @param fileB         File of a Matrix.
 */
static int read_column_ptr(uint64_t aCols, uint64_t** aColPtr, uint64_t bCols,
        uint64_t** bColPtr, FILE* fileA, FILE* fileB) {
    char* lineA = NULL;
    char* lineB = NULL;
    size_t lenA = 0;
    size_t lenB = 0;

    getline(&lineA, &lenA, fileA);
    getline(&lineB, &lenB, fileB);

    uint64_t pointerCountA = count_values(lineA);
    uint64_t pointerCountB = count_values(lineB);

    if (lineA && lineB) {
        if (pointerCountA != aCols + 1 || pointerCountB != bCols + 1) 
        {
            errno = EINVAL;
            perror("Wrong format. The number of column pointers is wrong");
            fprintf(stdout, "MatrixA Expected: %"PRIu64 " but was: %"PRIu64 "\n", pointerCountA, aCols + 1);
            fprintf(stdout, "MatrixB Expected: %"PRIu64 " but was: %"PRIu64 "\n", pointerCountB, bCols + 1);
            return 0;
        }
        uint64_t colPtrSizeA = aCols + 1;
        uint64_t colPtrSizeB = bCols + 1;

        *aColPtr = calloc(colPtrSizeA, sizeof(uint64_t));
        *bColPtr = calloc(colPtrSizeB, sizeof(uint64_t));
        if (!*aColPtr || !*bColPtr) {
            free(lineA);
            free(lineB);
            errno = ENOMEM;
			return 0;
        }

        // Store column Pointers
        if(line_parsing_uint(*aColPtr, colPtrSizeA, lineA, "column", 0) 
        || line_parsing_uint(*bColPtr, colPtrSizeB, lineB, "column", 0)) {
            free(lineA);
            free(lineB);
            return 1;
        }else {
            free(lineA);
            free(lineB);
            return 0;
        }
    } else {
        perror("An error occurred during parsing");
        free(lineA);
        free(lineB);
        return 0;
    }
}


int parse_csc_file(const char* filename_a, const char* filename_b, 
        struct cscMatrixTranspose* matrixA, struct cscMatrix* matrixB) {
    errno = 0;
    FILE* fileA = fopen(filename_a, "r");
    FILE* fileB = fopen(filename_b, "r");

    if (fileA == NULL || fileB == NULL) {
        perror("Unable to open file.");
        errno = ENOENT;
        return 0;    
    }

    // Read dimensions
    errno = 0;
    read_dimensions(&matrixA->rows, &matrixA->columns, &matrixB->rows, 
            &matrixB->columns, fileA, fileB);
    if (errno) return 0;

    // Read values
    read_values(&matrixA->valueCount, &matrixB->valueCount, &matrixA->values, 
            &matrixB->values, fileA, fileB);
    if (errno) return 0;

    // Read row indices
    read_row_indices(matrixA->valueCount, matrixA->rows, &matrixA->values, 
            &matrixA->rowIndices, matrixB->valueCount, matrixB->rows, &matrixB->values, 
            &matrixB->rowIndices, fileA, fileB);
    if (errno) return 0;

    // Read column pointers
    uint64_t isZero = read_column_ptr(matrixA->columns, &matrixA->colPtr, 
            matrixB->columns, &matrixB->colPtr, fileA, fileB);
    if (errno) {
        free(matrixA);
        free(matrixB);
        return 0;
    }
    
    //Get column indices 
    matrixA->colIndices = malloc(matrixA->valueCount * sizeof(uint64_t));
    if (matrixA->colIndices == NULL) {
        free_csc_matrix(matrixB);
        free_transpose_matrix(matrixA);
        errno = ENOMEM;
		return 0;
    }
    get_col_indices(matrixA, matrixA->colIndices);

    fclose(fileA);
    fclose(fileB);
    return isZero;
}


int parse_csc_file_V2(const char* filename_a, const char* filename_b, 
        struct cscMatrix* matrixA, struct cscMatrix* matrixB) {
    errno = 0;
    FILE* fileA = fopen(filename_a, "r");
    FILE* fileB = fopen(filename_b, "r");

    if (fileA == NULL || fileB == NULL) {
        perror("Unable to open file.");
        errno = ENOENT;
        return 0;    
    }

    // Read dimensions
    errno = 0;
    read_dimensions(&matrixA->rows, &matrixA->columns, &matrixB->rows, 
            &matrixB->columns, fileA, fileB);
    if (errno) return 0;

    // Read values
    read_values(&matrixA->valueCount, &matrixB->valueCount, &matrixA->values, 
            &matrixB->values, fileA, fileB);
    if (errno) return 0;

    // Read row indices
    read_row_indices(matrixA->valueCount, matrixA->rows, &matrixA->values, 
            &matrixA->rowIndices, matrixB->valueCount, matrixB->rows, &matrixB->values, 
            &matrixB->rowIndices, fileA, fileB);
    if (errno) return 0;

    // FIXME: (URGENT): matrix pointers used to be freed after some failures 
    // of read_column_ptr

    // Read column pointers
    uint64_t isZero = read_column_ptr(matrixA->columns, &matrixA->colPtr, 
            matrixB->columns, &matrixB->colPtr, fileA, fileB);
    if (errno) { // WARNING: (@Alejandro) I don't know if this is the right fix
        free(matrixA);
        free(matrixB);
        return 0;
    }

    fclose(fileA);
    fclose(fileB);
    return isZero;
}

void print_empty_matrix(uint64_t rows, uint64_t columns, const char* output_file) {
    FILE* file = fopen(output_file, "w");
    if (!file) {
        perror("Unable to open file ");
        errno = ENOENT;
		return;
    }

    //Write Dimensions
    fprintf(file,  "%"PRIu64 ",%"PRIu64 "\n", rows, columns);

    fprintf(file, "\n");
    fprintf(file, "\n");

    //Writte
    fprintf(file, "%u", 0);
    for (uint64_t i = 1; i < columns + 1; i++) {
        fprintf(file, ",%u", 0);
    }
    fclose(file);
}


void result_to_file(struct cscMatrix* result_matrix, const char* output_file) {  
    FILE* file = fopen(output_file, "w");
    if (!file) {
        perror("Unable to open file ");
        errno = ENOENT;
		return;
    }

    //Write Dimensions
    fprintf(file,  "%"PRIu64 ",%"PRIu64 "\n", result_matrix->rows, result_matrix->columns);
    //Write Values
    if (result_matrix->valueCount > 0) {
        fprintf(file, "%g", result_matrix->values[0]);

        for (uint64_t  i = 1; i < result_matrix->valueCount; i++) {   
            if (logData && !(i % 100)) {
                printf("\rWriting values. %.0f%% done.            ", 
                        100*((double) i)/result_matrix->valueCount);
                fflush(stdout);
            }

            fprintf(file, ",%g", result_matrix->values[i]);
        }
    }

    fprintf(file, "\n");

    if (result_matrix->valueCount > 0) {
        //Write row indices
        fprintf(file, "%"PRIu64 , result_matrix->rowIndices[0]);
        for (uint64_t i = 1; i < result_matrix->valueCount; i++) {
            if (logData && !(i % 100)) {
                printf("\rWriting row indices. "
                        "%.0f%% done.", 100*((double) i)/result_matrix->valueCount);
                fflush(stdout);
            }

            fprintf(file, ",%"PRIu64, result_matrix->rowIndices[i]);
        }
    }
    fprintf(file, "\n");

    //Write column pointers
    fprintf(file, "%" PRIu64, result_matrix->colPtr[0]);
    for (uint64_t i = 1; i < result_matrix->columns + 1; i++) {
            if (logData && !(i % 100)) {
            printf("\rWriting column pointers. "
                    "%.0f%% done.", 100*((double) i)/result_matrix->columns);
            fflush(stdout);
        }

        fprintf(file, ",%" PRIu64, result_matrix->colPtr[i]);
    }
    fclose(file);
}

