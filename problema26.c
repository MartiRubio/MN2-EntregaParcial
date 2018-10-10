#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/**
* Funció per a omplir el vector solució
*/
double* fill_vector(int n)
{
    double* vectorb = (double*)malloc(sizeof(double) * n);
    int i = 0;
    for(i = 0; i<n; i++){
        if (i%2 == 0){
            vectorb[i] = (i+2.)/n;
            vectorb[i+1] = (i+2.)/n;
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

/**
* Funció que retorna el valor de la matriu en la posició (i,j)
*/
int matrix_postion(int i, int j, int n)
{
    if (i == j){
        if(i%2==0){
            return 3;
        }
        else{
            return 4;
        }
    }
    else{
        if (abs(i-j) == 2){
            return -1;
        }
        else{
            if(abs(i-j) == n-1){
                return 1;
            }
            else{
                return 0;
            }
        }
    }
}

/**
* Funció que retorna el valor de la matriu en la posició (i,j)
*/
double matrix_bj_postion(int i, int j, int n)
{
    if(abs(i-j) == 2){
        if(i%2 == 0){
            return -1./3.;
        }
        else{
            return -1./4.;
        }
    }
    else{
        if(abs(i-j) == n-1){
            if(i%2 == 0){
                return 1./3.;
            }
            else{
                return 1./4.;
            }
        }
        else{
            return 0;
        }
    }
}

/**
* Funció que retorna el valor de la matriu en la posició (i,j)
*/
double vector_position(int i, int n)
{
    if(i%2 == 0){
        return (i+2.)/n;
    }
    else{
        return (i+1.)/n;
    }
}

int U_mes_L(int i, int j, int n){
    if (abs(i-j) == 2){
        return -1.;
    }
    else{
        if(abs(i-j) == n-1){
            return 1.;
        }
        else{
            return 0.;
        }
    }
}

double* D_inversa_mult(double* v, int n)
{
    double* vector;
    vector = fill_vector_solution(n);
    for(int i = 0; i < n; i++){
        if(i % 2 == 0){
            vector[i] = v[i]/3.;
        }
        else{
            vector[i] = v[i]/4.;
        }
    }
    return vector;
}

double* U_mes_L_mult(double* v, int n)
{
    double* vector;
    vector = fill_vector_solution(n);
    for(int i = 0; i < n; i++){
        if (i < 2){
            vector[i] += U_mes_L(i,i+2,n)*v[i+2];
            vector[i] += U_mes_L(i,n+1-i,n)*v[n+1-i];
        }
        else{
            if (i == n - 2){
                vector[i] += U_mes_L(i,i-2,n)*v[i-2];
                vector[i] += U_mes_L(i,0,n)*v[0];
            }
            else if (i == n-1){
                vector[i] += U_mes_L(i,i-2,n)*v[i-2];
                vector[i] += U_mes_L(i,1,n)*v[1];
            }
            else{
                vector[i] += U_mes_L(i,i-2,n)*v[i-2];
                vector[i] += U_mes_L(i,i+2,n)*v[i+2];
            }
        }
        
    }
    return vector;
}


void jacobi_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error, double modul_Bj)
{
    double error = 1;
    int Iteration = 1;
    printf("Iteration: %i\n", Iteration);
    // printf("%f\n", error);
    while (error > max_error){
        vector_solution_ant = vector_solution;
        vector_solution = U_mes_L_mult(vector_solution, n);
        for(int i = 0; i < n; i++){
            vector_solution[i] = vectorb[i] - vector_solution[i];
        }
        vector_solution = D_inversa_mult(vector_solution, n);
        error = infinite_norm(vector_solution, vector_solution_ant, n);
        // printf("E:    %f *10^-9\n", error*1000000000);
        Iteration++;
        printf("%i\n", Iteration);
    }
}


int main()
{
    // Creem les variables a usar
    int dimension = 1000000;
    double max_error = pow(10., -12.);
    double* vectorb;
    double* vector_solution;
    double* vector_solution_ant;
    double modul_Bj = 2./3.;

    // Omplim la matriu
    //matrixA = fill_matrix(dimension);

    // Omplim el vector
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);

    printf("\n\n********************\n");
    printf("*                  *\n");
    printf("* MÈTODE DE JACOBI *\n");
    printf("*                  *\n");
    printf("********************\n");
    jacobi_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error, modul_Bj);
    // Alliberem memòria

    printf("\n\n**************************\n");
    printf("*                        *\n");
    printf("* MÈTODE DE GAUSS-SEIDEL *\n");
    printf("*                        *\n");
    printf("**************************\n");
    /*
    for(i = 0; i < dimension; i++){
        free(matrixA[i]);
    }
    free(matrixA);
    */
    free(vectorb);
    free(vector_solution);
}
