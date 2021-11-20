#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

uint32_t binarySearch(uint32_t k, uint32_t i, uint32_t t_j, uint32_t const * const csc_col,uint32_t const * const csc_row);

void arrayMerging(uint32_t *array, uint32_t *subarr1, uint32_t *subarr2, uint32_t lensubarr1, uint32_t lensubarr2);

uint32_t commons(uint32_t *array1, uint32_t *array2, uint32_t lenarr1, uint32_t lenarr2);

#endif