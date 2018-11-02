#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define N 5

void jacobi(double *a, double *b, double *x, int n, int M) {
    double d, *u;
    u = (double*)calloc(n, sizeof(double));
    if(u == NULL) {
        fprintf(stderr, "Error while allocating memory\n");
        exit(1);
    }
    for(int i = 0; i < n; i++) {
        d = 1 / a[i*n+i];
        b[i] *= d;
        for(int j = 0; j < n; j++) {
            a[i*n+j] *= d;
        }
    }
    for(int k = 0; k < M; k++) {
        for(int i = 0; i < n; i++) {
            double sum = 0.0;
            for(int j = 0; j < n; j++) {
                if(j == i) 
                    continue;
                else 
                    sum += a[i*n+j] * x[j];
            }
            u[i] = b[i] - sum;
        }
        for(int i = 0; i < n; i++) {
            x[i] = u[i];
        }
        printf("%d. Jacobi iteration result:\n", k);
        for(int i = 0; i < n; i++) {
            printf("%f ", x[i]);
        }
        printf("\n\n");
    }
    free(u);
}

void gaussian_elimination(double *a, double *b, double *x, int n) {

    int column, row, diagonal, max_pivot_row, j;
    double max_pivot, tmp;
    for(diagonal = 0; diagonal < n; diagonal++) {
        max_pivot_row = diagonal;
        max_pivot = *(a + (diagonal * n + diagonal));    // i,ith element of the matrix
        for(row = diagonal + 1; row < n; row++) {
            tmp = fabs(*(a + (row * n + diagonal)));
            if(tmp > max_pivot) {
                max_pivot_row = row;
                max_pivot = tmp;
            }
        }

        if(diagonal != max_pivot_row) {
            for(int k = 0; k < n; k++) {
                double *tmp_pointer1 = a + (diagonal * n + k);
                double *tmp_pointer2 = a + (max_pivot_row * n + k);
                tmp = *tmp_pointer1;
                *tmp_pointer1 = *tmp_pointer2;
                *tmp_pointer2 = tmp;
            }
            tmp = b[diagonal];
            b[diagonal] = b[max_pivot_row];
            b[max_pivot_row] = tmp;
        }

        for(row = diagonal + 1; row < n; row++) {
            tmp = *(a + (row * n + diagonal)) / *(a + (diagonal * n + diagonal));
            for(column = diagonal + 1; column < n; column++) {
                *(a + (row * n + column)) -= tmp * *(a + (diagonal * n + column));
            }
            *(a + (row * n + diagonal)) = 0;
            b[row] -= tmp * b[diagonal];
        }
    }

    for(row = n - 1; row >= 0; row--) {
        tmp = b[row];
        for(j = n - 1; j > row; j--) {
            tmp -= x[j] * *(a + (row * n + j));
        }
        x[row] = tmp / *(a + (row * n + row));
    }

}

int main(void) {

    static double a[N * N], a_prim[N * N];
    static double b[N], b_prim[N];
    static double x[N], x_prim[N];
    int s;

    FILE *a_file = fopen("1DLaplacianA.txt", "r");
    FILE *b_file = fopen("1DLaplacianB.txt", "r");
    if(a_file == NULL || b_file == NULL) {
        fprintf(stderr, "Error while opening files.\n");
        exit(1);
    }

    srand(time(NULL));

    printf("Matrix A:\n");
    for(int i = 0; i < N * N; i++) {
#ifndef RND
        fscanf(a_file ,"%lf", &a[i]);
#else
        a[i] = (double)random()/RAND_MAX * 5;
#endif
        a_prim[i] = a[i];
        printf("%g ", a[i]);
        if((i+1) % N == 0 && i) {
            printf("\n");
        }
    }
    printf("\n");

    printf("Matrix B:\n");
    for(int i = 0; i < N; i++) {
#ifndef RND
        fscanf(b_file ,"%lf", &b[i]);
#else
        b[i] = (double)random()/RAND_MAX * 5;
#endif
        b_prim[i] = b[i];
        x[i] = 0.0;
        x_prim[i] = 0.0;
        printf("%g\n", b[i]);
    }
    printf("\n");

    jacobi(a, b, x, N, 70);

    gaussian_elimination(a_prim, b_prim, x_prim, N);

    printf("\nGaussian elimination with pivoting result:\n");
    for(int i = 0; i < N; i++) {
        printf("%.15g\n", x_prim[i]);
    }
    fclose(a_file);
    fclose(b_file);

    return 0;

}