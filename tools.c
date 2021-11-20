#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void arrayMerging(uint32_t *array, uint32_t *subarray1, uint32_t *subarray2, uint32_t lensubarr1, uint32_t lensubarr2) {
    uint32_t i = 0;
    uint32_t j = 0;
    while (i < lensubarr1 && j < lensubarr2) {
        if (subarray1[i] < subarray2[j]) {
            array[i + j] = subarray1[i];
            i++;
        } else {
            array[i + j] = subarray2[j];
            j++;
        }
    }
    while (i < lensubarr1) {
        array[i + j] = subarray1[i];
        i++;
    }
    while (j < lensubarr2) {
        array[i + j] = subarray2[j];
        j++;
    }
}

uint32_t commons(uint32_t *array1, uint32_t *array2, uint32_t lenarr1, uint32_t lenarr2) {
    uint32_t total = 0;
    uint32_t i = 0;
    uint32_t j = 0;
    while (i < lenarr1 && j < lenarr2) {
        if (array1[i] < array2[j]) {
            i++;
        } else if (array1[i] > array2[j]) {
            j++;
        } else {
            i++;
            j++;
            total++;
        }
    }
    return total;
}