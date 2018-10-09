#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/**
* Funció per a omplir la matriu
*/
void fill_matrix(double **matrix, int n)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < n; i++) {
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
    fill_diagonal(matrix, n);
    return matrix;
}

/**
* Funció per a omplir la diagonal de la matriu
*/
void fill_diagonal(double **matrix, int n)
{
    int i = 0;
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
void fill_vector(double **matrix)
{
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
* Funció per a calcular la norma infinit d'una matriu
*/
void infinite_norm(double **matrix)
{

}

int main()
{
    // Creem les variables a usar
    int dimension = 1000000;
    double max_error = pow(10., -12.);
    double* vectorb = (double*)malloc(sizeof(double) * dimension);
    double* vector_solution = (double*)malloc(sizeof(double) * dimension);
    double** matrixA = (double**)malloc(sizeof(double*) * dimension);
    
    // Omplim la matriu
    fill_matrix(matrixA, dimension);

    // Omplim el vector
    fill_vector(vectorb, dimension);
    

    // Alliberem memòria
    int i = 0;
    for(i = 0; i < dimension; i++){
        free(matrixA[i]);
    }
    free(matrixA);
    free(vectorb);
    free(vector_solution);
}
