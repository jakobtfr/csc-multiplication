#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "csc_io.h"
#include "matrix_mul.h"
#include "cs_matrix.h"
#include "transpose.h"

static double get_time_diff(struct timespec* start, struct timespec* end) {
    return end->tv_sec - start->tv_sec + 
        1e-9 * (end->tv_nsec - start->tv_nsec);
}

static void get_time(struct timespec* t) {
    clock_gettime(CLOCK_MONOTONIC, t);
}

int main(int argc, char *argv[]) {
    const char* progname = argv[0];

    int version = 0;
    char* file_a = "data/matrixA.txt";
    char* file_b = "data/matrixB.txt";
    char* output_file = "data/result.txt";
    unsigned int iterations = 1;
    int measureTime = 0;
    logData = 0;
    int option_index = 0;
    int generateNew = 0;

    int opt;
    while ((opt = getopt_long(argc, argv, shortopts, longopts, &option_index)) != -1) {
        switch (opt) {
            case 'V':
                if (convert_int(optarg, &version) != 0) {
                    return EXIT_FAILURE;
                }
                break;
            case 'a':
                file_a = optarg;
                break;
            case 'b':
                file_b = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'B':
                measureTime = 1;
                if (optarg) {
                    
                    if (convert_unsigned(optarg, &iterations) != 0) {
                        return EXIT_FAILURE;
                    }
                }
                break;
            case 'h':
                print_help(progname);
                return EXIT_SUCCESS;
            case 'l':
                logData = 1;
                break;
            case 'r':
                generateNew = 1;
                break;
            default:
                abort();
        }
    }

    srand(time(NULL));

    int (*transpose_fun)(struct cscMatrixTranspose*, struct cscMatrix*);
    void (*mul_fun)(const void*, const void*, void*);

    switch (version) {
        case 0:
            transpose_fun = transpose;
            mul_fun = matr_mult_csc;
            break;
        case 1:
            transpose_fun = transpose;
            mul_fun = matr_mult_csc_V1;
            break;
        case 2:
            transpose_fun = 0;
            mul_fun = matr_mult_csc_V2;
            break;
        default:
            fprintf(stderr, "Invalid program version.\n");
            return EXIT_FAILURE;
            break;
    }

    struct timespec start_time, end_time, mul_start, mul_end, 
            transpose_start, transpose_end, create_start, create_end, 
            parse_start, parse_end;
    double elapsed_time, mul_time = 0, transpose_time = 0, create_time, 
          parse_time = 0;
    if (measureTime) get_time(&start_time);

    if (generateNew) {
        if (measureTime) get_time(&create_start);
        errno = 0;
        create_random_matrix_doc(512, 1024, "randomMatrixA.txt", "randomMatrixB.txt");

        int err = errno;
        if (measureTime) get_time(&create_end);
        if (err) return 1;
    }

    for (size_t i = 0; i < iterations; i++) {

        struct cscMatrix* matrixB = malloc(sizeof(struct cscMatrix));
        struct cscMatrix* matrixA_in = malloc(sizeof(struct cscMatrix));
        if (!matrixA_in || !matrixB) {
            errno = ENOMEM;
            perror("Error initializing input matrices");
            return EXIT_FAILURE;
        }
        struct cscMatrix* result = malloc(sizeof(struct cscMatrix));
        if (!result) {
            perror("Error initializing output matrix");
            return EXIT_FAILURE;
        }

        if (logData) printf("Input matrix structs initialized successfully.\n");

        if (measureTime) get_time(&parse_start);
        errno = 0;
        if (version != 2) {
            struct cscMatrixTranspose* matrixA = 
                malloc(sizeof(struct cscMatrixTranspose));
            if (!matrixA) {
                perror("Error initializing auxiliary matrix");
                return EXIT_FAILURE;
            }
            if(parse_csc_file(file_a, file_b, matrixA, matrixB)){
                if (measureTime) {
                    get_time(&parse_end);
                    parse_time += get_time_diff(&parse_start, &parse_end);
                }
                print_empty_matrix(matrixA->columns, matrixB->columns, 
                        output_file);
                free_csc_matrix(matrixB);
                free_transpose_matrix(matrixA);
                continue;
            }
            if (measureTime) {
                get_time(&parse_end);
                parse_time += get_time_diff(&parse_start, &parse_end);
            }


            if (errno || !matrixA || !matrixB) {
                fprintf(stderr, "Matrix parsing failed.\n");
                return EXIT_FAILURE;
            }



            if (logData) printf("Transposing matrix A...");
            if (measureTime) get_time(&transpose_start);
            if (!transpose_fun(matrixA, matrixA_in)) {
                if (measureTime) get_time(&transpose_end);
                if (logData) perror("\rError transposing matrix A");
                free_transpose_matrix(matrixA);
                free_csc_matrix(matrixB);
                return EXIT_FAILURE;
            } else {
                free_transpose_matrix(matrixA);
            }
        } else {
            if(parse_csc_file_V2(file_a, file_b, matrixA_in, matrixB)){
                if (measureTime) {
                    get_time(&parse_end);
                    parse_time += get_time_diff(&parse_start, &parse_end);
                }
                print_empty_matrix(matrixA_in->rows, matrixB->columns,
                        output_file);
                free_csc_matrix(matrixB);
                free_csc_matrix(matrixA_in);
                free_csc_matrix(result);
                continue;
            }
            if (measureTime) {
                get_time(&parse_end);
                parse_time += get_time_diff(&parse_start, &parse_end);
            }


            if (errno || !matrixA_in || !matrixB) {
                fprintf(stderr, "Matrix parsing failed.\n");
                return EXIT_FAILURE;
            }
        }

        if (measureTime) {
            get_time(&transpose_end);
            transpose_time += get_time_diff(&transpose_start, &transpose_end);
        }

        if(logData) printf("\rMatrix A transposed successfully.\n");

        if (measureTime) get_time(&mul_start);
        mul_fun(matrixA_in, matrixB, result); 
        // Store dimensions and valueCounts for logging
        uint64_t aRows = matrixA_in->rows, aCols = matrixA_in->columns,
                 aVals = matrixA_in->valueCount, bRows = matrixB->rows,
                 bCols = matrixB->columns, bVals = matrixB->valueCount;
        free_csc_matrix(matrixB);
        free_csc_matrix(matrixA_in);
        if (measureTime) {
            get_time(&mul_end);
            mul_time += get_time_diff(&mul_start, &mul_end);
        }

        if (errno != 0) {
            fprintf(stderr, "Matrix multiplication failed.\n");
            free(result);
            return EXIT_FAILURE;
        }
        
        
        if (logData) {
            printf("Writing result to output file...");
            fflush(stdout);
        }
        if (measureTime) {
            get_time(&parse_start);
        }
        result_to_file(result, output_file);
        if (measureTime) {
            get_time(&parse_end);
            parse_time += get_time_diff(&parse_start, &parse_end);
        }

        if (errno != 0) {
            perror("\rError writing result to output file");
            return EXIT_FAILURE;
        }
        if (logData && i == iterations - 1) {
            printf("\rOutput written to file successfully.\n"
                "Operation %lu complete.\n", i+1);
            printf("Data logs:\n");
            printf("Matrix A\n");
            printf("Dimensions: %lu rows and %lu columns\n", 
                    aCols, aRows);
            printf("Density: %.1f%%\n\n", 
                    100*((double) aVals)/
                    (aCols * aRows));
            printf("Matrix B\n");
            printf("Dimensions: %lu rows and %lu columns\n", 
                    bRows, bCols);
            printf("Density: %.1f%%\n\n", 
                    100*((double) bVals)/
                    (bCols * bRows));
            printf("Result Matrix\n");
            printf("Dimensions: %lu rows and %lu columns\n", 
                    result->rows, result->columns);
            printf("Density: %.1f%%\n\n", 
                    100*((double) result->valueCount)/
                    (result->rows*result->columns));
        }

        // Clean up
        free_csc_matrix(result);
        
    }

    if (iterations == 1){
        printf("The program has been run %u time.\n", iterations);
    }else {
        printf("The program has been run %u times.\n", iterations);
    }
    printf("Version: %d\n", version);

    if (measureTime)
    {
        get_time(&end_time); // Stop the timer
        elapsed_time = get_time_diff(&start_time, &end_time);
        if (generateNew) {
            create_time = get_time_diff(&create_start, &create_end);
        } else create_time = 0;
        double average_create = create_time / iterations;

        if (iterations > 1) {
            printf("Average computation time: %g s.\n", (mul_time / iterations));
            if (version != 2) {
                printf("Average transposing time: %g s.\n", 
                        (transpose_time / iterations));
            }
            printf("Average parsing time: %g s.\n", parse_time / iterations);
            printf("Average iteration length: %g s.\n", 
                    (elapsed_time / iterations) - average_create);
            printf("Total time for %u iterations: %g seconds.\n", iterations, 
                    elapsed_time - create_time);
        } else {
            printf("Total time: %g seconds.\n", elapsed_time);
        }

        printf("Total I/O processing time: %g s.\n", parse_time);
        printf("Total computation time: %g s.\n", mul_time);
        if (version != 2) {
            printf("Total transposing time: %g s.\n", transpose_time);
        }
    }
    return EXIT_SUCCESS;
}

