#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


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
        vectorx[i] = 0.;
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
            return 3.;
        }
        else{
            return 4.;
        }
    }
    else{
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
* Funció que retorna el valor del vector b en la fila i
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

/*
* Funció que multiplica el vector v per la matriu D^{-1}
*/
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

/*
* Funció que multiplica el vector v per la matriu U+L
*/
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

double U_mes_L_fila_i_mult(double* v, double* v_ant, int i, int n){
    double result;
    if(i < 2){
        result = v_ant[i + 2]*U_mes_L(i, i+2, n) + v_ant[n-2+i]*U_mes_L(i, n-2+i, n);
    }
    else if(i > n - 2){
        result = v[i - 2]*U_mes_L(i, i-2, n) + v[i-n+2]*U_mes_L(i, i-n+2, n);
    }
    else{
        result = v[i - 2]*U_mes_L(i, i-2, n) + v_ant[i + 2]*U_mes_L(i, i+2, n);
    }
    return result;
}


void assign_vectors(double* v_ant, double* v, int n){
    for(int i = 0; i < n; i++){
        v_ant[i] = v[i];
    }
}

/*
* Mètode que aplica l'algoritme de Jacobi
*/
void jacobi_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error)
{
    double error = 1.;
    int Iteration = 1;
    // printf("Iteration: %i\n", Iteration);
    // printf("%f\n", error);
    while (error > max_error){
        //assign_vectors(vector_solution_ant, vector_solution, n);
        vector_solution_ant = vector_solution;
        vector_solution = U_mes_L_mult(vector_solution, n);
        for(int i = 0; i < n; i++){
            vector_solution[i] = vectorb[i] - vector_solution[i];
        }
        vector_solution = D_inversa_mult(vector_solution, n);
        error = infinite_norm(vector_solution, vector_solution_ant, n);
        // printf("E:    %f *10^-9\n", error*1000000000);
        Iteration++;
        // printf("%i\n", Iteration);
    }
    printf("Total d'iteracions: %i\n", Iteration - 1);
}

/*
* Mètode que aplica l'algoritme de Gauss-Seidel
*/
void gauss_seidel_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error){
    double error = 1.;
    int Iteration = 1;
    // printf("Iteration: %i\n", Iteration);
    // printf("%f\n", error);
    while (error > max_error){
        assign_vectors(vector_solution_ant, vector_solution, n);
        // vector_solution_ant = vector_solution;
        for(int i = 0; i < n; i++){
            vector_solution[i] = vectorb[i] - U_mes_L_fila_i_mult(vector_solution, vector_solution_ant, i, n);
            vector_solution[i] = vector_solution[i]/matrix_postion(i,i,n);
        }
        error = infinite_norm(vector_solution, vector_solution_ant, n);
        // printf("E:    %f *10^-9\n", error*1000000000);
        Iteration++;
        // printf("%i\n", Iteration);
    }
    printf("Total d'iteracions: %i\n", Iteration - 1);
}

/*
* Mètode que aplica l'algoritme de Gauss-Seidel
*/
void SOR_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error, double omega){
    double error = 1.;
    int Iteration = 1;
    printf("Iteration: %i\n", Iteration);
    printf("%f\n", error);
    while (error > max_error){
        assign_vectors(vector_solution_ant, vector_solution, n);
        // vector_solution_ant = vector_solution;
        for(int i = 0; i < n; i++){
            vector_solution[i] = vectorb[i] - U_mes_L_fila_i_mult(vector_solution, vector_solution_ant, i, n);
            vector_solution[i] = vector_solution[i]/matrix_postion(i,i,n);
            vector_solution[i] = vector_solution_ant[i] + vector_solution[i]*omega;
        }
        error = infinite_norm(vector_solution, vector_solution_ant, n);
        printf("E:    %f *10^-9\n", error*1000000000);
        Iteration++;
        printf("%i\n", Iteration);
    }
    printf("Total d'iteracions: %i\n", Iteration - 1);
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
    double omega = -0.2;

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
    
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    jacobi_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode de Jacobi: %f\n", cpu_time_used);

    // Tornem a omplir els vectors 
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);
    printf("\n\n**************************\n");
    printf("*                        *\n");
    printf("* MÈTODE DE GAUSS-SEIDEL *\n");
    printf("*                        *\n");
    printf("**************************\n");

    start = clock();
    gauss_seidel_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode de Gauss-Seidel: %f\n", cpu_time_used);

    // Tornem a omplir els vectors 
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);
    printf("\n\n**************\n");
    printf("*            *\n");
    printf("* MÈTODE SOR *\n");
    printf("*            *\n");
    printf("**************\n");

    start = clock();
    SOR_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error, omega);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode SOR: %f\n", cpu_time_used);
    /*
    for(i = 0; i < dimension; i++){
        free(matrixA[i]);
    }
    free(matrixA);
    */
    free(vectorb);
    free(vector_solution);
}
