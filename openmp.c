#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "mmio.h"
#include "coo2csc.h"

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

int main(int argc, char *argv[]) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nnz;   
    int *I, *J;
    double *val;
    int threads_number = atoi(argv[2]);
    char* threads_number_char = argv[2];
    struct timeval start, end;

    if (argc < 2) {
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	} else {
        if ((f = fopen(argv[1], "r")) == NULL) {
            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nnz)) != 0) {
        exit(1);
    }


    // Reserving memory for COO, CSC, C_CSC matrices
    I = (uint32_t *) malloc(2 * nnz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2 * nnz * sizeof(uint32_t));
    val = (double *) malloc(nnz * sizeof(double));
    uint32_t* csc_row = (uint32_t *) malloc(2 * nnz * sizeof(uint32_t));
    uint32_t* csc_col = (uint32_t *) malloc((N + 1) * sizeof(uint32_t));
    uint32_t* c_csc_row = (uint32_t *) malloc(0 * sizeof(uint32_t));
    uint32_t* c_val = (uint32_t *) malloc(0 * sizeof(uint32_t));
    uint32_t* c_csc_col = (uint32_t *) malloc((N + 1) * sizeof(uint32_t));

    for (uint32_t i = 0; i < nnz; i++) {
            /* I is for the rows and J for the columns */
            fscanf(f, "%d %d \n", &I[i], &J[i]);
            I[i]--;  /* adjust from 1-based to 0-based */
            J[i]--;
    }
    
    if (f != stdin) {
        fclose(f);
    }

    if (M != N) {
        printf("Columns and rows differ in size.");
    }

    // Generating symmetrical matrix
    for (uint32_t i = 0; i < nnz; i++) {
        I[nnz + i] = J[i];
        J[nnz + i] = I[i];
    }

    // Swapping I and J according to the symmetrical matrix, to achieve an upper triangular matrix
    if (I[0] > J[0]) {
        coo2csc(csc_row, csc_col, J, I, 2 * nnz, M, 0);
    } else {
        coo2csc(csc_row, csc_col, I, J, 2 * nnz, N, 0);
    }

    printf("\nMatrix Loaded!\n");

    // Initializing c3 and results with zeros and e with ones
    int *c3, *e, *results;
    c3 = malloc(N * sizeof c3);    
    e = malloc(N * sizeof e);  
    results = malloc(N * sizeof results);  
    for (int i = 0; i < N; i++){
        c3[i] = 0;
        e[i] = 1;
        results[i] = 0;
    }

    c_csc_col[0] = 0;
    c_csc_row = realloc(c_csc_row, 2 * nnz * sizeof(int));
    c_val = realloc(c_csc_row, 2 * nnz * sizeof(int)); 
    
    // Time measurements starts here, as we begin the sequence of calculating the Hadamard product and the triangles
    gettimeofday(&start,NULL);  
    int *l;

    // Hadamard Product. We use threads_number for all regions, dynamic teams disabled
    omp_set_dynamic(0);
    omp_set_num_threads(threads_number);
    
    #pragma omp parallel for private(l)
    for (int i = 0; i < N; i++) {
        for(int j = 0; j < csc_col[i+1] - csc_col[i]; j++) {
            int current_row = csc_row[csc_col[i] + j];
            int current_col = i;

            // Calculating A*A
            int alpha_size = csc_col[current_row+1] - csc_col[current_row];
            int *alpha = malloc((alpha_size) * sizeof(int));
            int beta_size = csc_col[current_col+1] - csc_col[current_col];    
            int *beta = malloc((beta_size) * sizeof(int));
            for(int k = 0; k < alpha_size; k++) {
                alpha[k] = csc_row[csc_col[current_row] + k];
            }
            for(int k = 0; k < beta_size; k++) {
                beta[k] = csc_row[csc_col[current_col] + k];
            }

            int value = commons(alpha, beta, alpha_size, beta_size);       

            if(value) {
                c_val[csc_col[i] + j] = value;
            }
            free(beta);
            free(alpha);
        }
    }
    c_csc_row = csc_row;
    c_csc_col = csc_col;

    //Final Multiplication - no profit taken out of parallelization
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < c_csc_col[i+1] - c_csc_col[i]; j++) {
            int row = c_csc_row[c_csc_col[i] + j];
            int col = i;
            int value = c_val[c_csc_col[i] + j];
            #pragma omp critical
            results[row] += value * e[col];
        }
    }
    int totalTriangles = 0;
    for(int i = 0; i < N; i++) {
        c3[i] = results[i] / 2;
        totalTriangles += c3[i];
    }

    totalTriangles = totalTriangles / 3;

    // End of procedure and 'end' timestamp
    gettimeofday(&end,NULL);
    double duration = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    printf("Result: %d triangles total.\n",  totalTriangles);
    printf("Duration: %f seconds.\n",  duration);

    // Free arrays
    free(I);
    free(J);
    free(c3);
    free(e);
    free(results);

	return 0;
}