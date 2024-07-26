#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "cs_matrix.h"
#include "radixsort.h"
/**
 * utility function to get max value in arr
 * @param arr   array to get max of
 * @param n     size of array
 */
static uint64_t getMax(uint64_t* arr, uint64_t n) {
    if (n < 1) {
        errno = EINVAL;
        return 0;
    }
    uint64_t max = *arr;
    for (uint64_t i = 1; i < n; ++i) {
        if (*(arr + i) > max)
            max = *(arr + i);
    }
    return max;
}
/**
 * Computes the next greater equal power of 2 of a given value
 * @param x     The given 64bit value
 * @return      Next power of 2 greater equal to x
 */
/*this method is currently not used, for eventual binary Radixsort
static uint64_t next_pow2(uint64_t x) {
    x |= x>>1;
    x |= x>>2;
    x |= x>>4;
    x |= x>>8;
    x |= x>>16;
    x |= x>>32;
    return x + 1;
    // https://graphics.stanford.edu/~seander/bithacks.html
    // extended method to be able to handle 64bit values
} */


/**
 * Sorts elements of arr based on significant places
 * https://de.wikipedia.org/wiki/Countingsort implementation of pseudocode
 * @param keys   Array to be sorted
 * @param val1   first Array that should be sorted by keys
 * @param val2   first Array that should be sorted by keys
 * @param n      size of keys
 * @param place  current significant digit
 */
static void kSort(uint64_t* keys, float* val1, uint64_t* val2, uint64_t n, uint64_t exp) {
    uint64_t* outKeys = malloc(n * sizeof(uint64_t));
    float* outVal1 = malloc(n * sizeof(float));
    uint64_t* outVal2 = malloc(n * sizeof(uint64_t));
    uint64_t count[10] = {0};

    if (!outKeys || !outVal1 || !outVal2) {
        errno = ENOMEM;
        free(outKeys);
        free(outVal1);
        free(outVal2);
        return;
    }

    // get element count
    for (uint64_t i = 0; i < n; ++i)
        count[(keys[i] / exp) % 10]++;


    for (uint64_t i = 1; i < 10; ++i)
        count[i] += count[i - 1];


    // sort elements
    for (uint64_t i = n - 1; i > 0; --i) {
        uint64_t index = (keys[i] / exp) % 10;
        outKeys[count[index] - 1] = keys[i];
        outVal1[count[index] - 1] = val1[i];
        outVal2[count[index] - 1] = val2[i];
        count[index]--;
    }
    uint64_t index = (keys[0] / exp) % 10;
    outKeys[count[index] - 1] = keys[0];
    outVal1[count[index] - 1] = val1[0];
    outVal2[count[index] - 1] = val2[0];
    count[index]--;
    // sort key array and sort val1 and val2 by that order
    for (uint64_t i = 0; i < n; i++) {
        keys[i] = outKeys[i];
        val1[i] = outVal1[i];
        val2[i] = outVal2[i];
    }
    free(outKeys);
    free(outVal1);
    free(outVal2);

}

/**
 * Main function that sorts arr with radixSort
 * @param keys      array to be sorted
 * @param values    first Array that should be sorted by keys
 * @param val2      second Array that should be sorted by keys
 * @param n         size of the array
 */
void radixSort(uint64_t* keys, float* values, uint64_t* val2, uint64_t size) {
    errno = 0;
    uint64_t max = getMax(keys, size);
    if (max == 0 && errno == EINVAL) return;
    for (uint64_t place = 1; max / place > 0; place *= 10) {
        kSort(keys,  values, val2, size, place);
        if (errno == ENOMEM) {
            perror("Could not allocate memory for transposing procedure");
        }
    }

}
