#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/**
* Funció per a omplir la matriu
*/
double** fill_matrix(int n)
{
    int i = 0;
    int j = 0;
    double** matrix = (double**)malloc(sizeof(double*) * n);
    for (i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(sizeof(double) * n);
        for (j = 0; j < n; j++) {
            if((i+j) % 2 != 0){
                matrix[i][j] = 0;
            }
            else{
                if((i+j) % 4 == 0){
                    matrix[i][j] = 1;
                }
                else{
                    matrix[i][j] = -1;  
                }
            }
        }
    }
    // Omplim la diagonal
    for(i=0; i < n; i++){
        if(i%2==0){
            matrix[i][i] = 3;
        }
        else{
            matrix[i][i] = 4;
        }
    }
    return matrix;
}

/**
* Funció per a omplir el vector
*/
double* fill_vector(int n)
{
    double* vectorb = (double*)malloc(sizeof(double) * n);
    int i = 0;
    for(i = 0; i<n; i++){
        if(i%2==0){
            vectorb[i] = i+2/n;
            vectorb[i+1] = i+2/n;
        }
    }
    return vectorb;
}

/**
* Funció per a omplir el vector solució
*/
double* fill_vector_solution(int n)
{
    double* vectorx = (double*)malloc(sizeof(double) * n);
    int i = 0;
    for(i = 0; i<n; i++){
        vectorx[i] = 0;
    }
    return vectorx;
}

/**
* Funció per a calcular la norma infinit d'una matriu
*/
double infinite_norm(double *vectorb, double *vectorx, int n)
{
    int i = 1;
    double norm = fabs(vectorb[0] - vectorx[0]);
    for (i = 1; i <n; i++){
        if (norm > fabs(vectorb[i] - vectorx[i])){
            norm = fabs(vectorb[i] - vectorx[i]);
        }
    }
    return norm;
}

int main()
{
    // Creem les variables a usar
    int dimension = 1000000;
    double max_error = pow(10., -12.);
    double* vectorb;
    double* vector_solution;
    double** matrixA;
    
    // Omplim la matriu
    matrixA = fill_matrix(dimension);

    // Omplim el vector
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);

    // Alliberem memòria
    int i = 0;
    for(i = 0; i < dimension; i++){
        free(matrixA[i]);
    }
    free(matrixA);
    free(vectorb);
    free(vector_solution);
}
