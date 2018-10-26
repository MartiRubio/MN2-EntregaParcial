// Author: Martí Rubio
// Codi @ https://github.com/MartiRubio/MN2-EntregaParcial

/*
Codi per al problema 26 de la llista 1 de Mètodes Numèrics II
*/

// Includes necessaris
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


//Functions declaration
double* fill_vector(int);
double* fill_vector_solution(int);
double infinite_norm(double*, double*, int);
int matrix_postion(int, int, int);
int U_mes_L(int, int, int);
double* D_inversa_mult(double*, int);
double* U_mes_L_mult(double*, int);
double U_mes_L_fila_i_mult(double*, double*, int, int);
double U_mes_L_fila_i_mult(double*, double*, int, int);
void assign_vectors(double*, double*, int);
void jacobi_method(double*, double*, double*, int, double);
void gauss_seidel_method(double*, double*, double*, int, double);
void SOR_method(double*, double*, double*, int, double, double);


/**
* Funció per a omplir el vector b
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
* Funció per a omplir el vector solució (inicialitzat a zero)
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
* Funció per a calcular la norma infinit entre dos vectors
*/
double infinite_norm(double *vectorb, double *vectorx, int n)
{
    int i = 1;
    double norm = fabs(vectorb[0] - vectorx[0]);
    for (i = 1; i <n; i++){
        if (norm < fabs(vectorb[i] - vectorx[i])){
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
            if(abs(i-j) == n-2){
                return 1.;
            }
            else{
                return 0.;
            }
        }
    }
}


/**
* Funció que retorna el valor de la matriu U+L a la posició (i,j)
*/
int U_mes_L(int i, int j, int n)
{
    if (abs(i-j) == 2){
        return -1.;
    }
    else{
        if(abs(i-j) == n-2){
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
            vector[i] += U_mes_L(i,n-2+i,n)*v[n-2+i];
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


/*
* Funció que multiplica el vector v i v_anterior per la matriu U+L (mode Gauss Seidel)
*/
double U_mes_L_fila_i_mult(double* v, double* v_ant, int i, int n)
{
    double result = 0.;
    if(i < 2){
        result = v_ant[i + 2]*U_mes_L(i, i+2, n) + v_ant[n-2+i]*U_mes_L(i, n-2+i, n);
    }
    else if(i > n - 3){
        result = v[i - 2]*U_mes_L(i, i-2, n) + v[i-n+2]*U_mes_L(i, i-n+2, n);
    }
    else{
        result = v[i - 2]*U_mes_L(i, i-2, n) + v_ant[i + 2]*U_mes_L(i, i+2, n);
    }
    return result;
}


/*
* Posem al vector anterior el contingut del vector actual
*/
void assign_vectors(double* v_ant, double* v, int n)
{
    for(int i = 0; i < n; i++){
        v_ant[i] = v[i];
    }
}

/*
* Mètode que aplica l'algoritme de Jacobi
*/
void jacobi_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error)
{

    // variables de la fució
    double error = 1.;
    int Iteration = 1;

    // Mentre l'error sigui major al mínim esperat
    while (error > max_error){

        // Posem al vector anterior el vector actual
        vector_solution_ant = vector_solution;

        // Apliquem el mètode de Jacobi
        vector_solution = U_mes_L_mult(vector_solution, n);
        for(int i = 0; i < n; i++){
            vector_solution[i] = vectorb[i] - vector_solution[i];
        }
        vector_solution = D_inversa_mult(vector_solution, n);

        // Calculem la norma infinit entre els dos vectors
        error = infinite_norm(vector_solution, vector_solution_ant, n);

        // Calculem la norma infinit entre els dos vectors
        Iteration++;
    }

    // Imprimim el nombre d'iteracions
    printf("Total d'iteracions: %i\n", Iteration - 1);
}


/*
* Mètode que aplica l'algoritme de Gauss-Seidel
*/
void gauss_seidel_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error)
{

    // Variables de la funció
    double error = 1.;
    int Iteration = 1;

    // Mentre l'error sigui major al mínim esperat
    while (error > max_error){

        // Posem al vector anterior el vector actual
        assign_vectors(vector_solution_ant, vector_solution, n);

        // Per cada fila de la matriu
        for(int i = 0; i < n; i++){

            // Apliquem el mètode de Gauss-Seidel
            vector_solution[i] = vectorb[i] - U_mes_L_fila_i_mult(vector_solution, vector_solution_ant, i, n);
            vector_solution[i] = vector_solution[i]/matrix_postion(i,i,n);
        }

        // Calculem la norma infinit entre els dos vectors
        error = infinite_norm(vector_solution, vector_solution_ant, n);

        // Sumem una iteració
        Iteration++;
    }

    // Imprimim el nombre d'iteracions
    printf("Total d'iteracions: %i\n", Iteration - 1);
}


/*
* Mètode que aplica l'algoritme SOR
*/
void SOR_method(double *vector_solution, double *vector_solution_ant, double *vectorb, int n, double max_error, double omega)
{

    // Variables de la funció
    double error = 1.;
    int Iteration = 1;

    // Mentre l'error sigui major al mínim esperat
    while (error > max_error){

        // Posem al vector anterior el vector actual
        assign_vectors(vector_solution_ant, vector_solution, n);

        // Per cada fila de la matriu
        for(int i = 0; i < n; i++){

            // Apliquem el mètode SOR
            vector_solution[i] = vectorb[i] - U_mes_L_fila_i_mult(vector_solution, vector_solution_ant, i, n);
            vector_solution[i] = vector_solution[i]/matrix_postion(i,i,n);
            vector_solution[i] = (1 - omega)*vector_solution_ant[i] + vector_solution[i]*omega;
        }

        // Calculem la norma infinit entre els dos vectors
        error = infinite_norm(vector_solution, vector_solution_ant, n);

        // Sumem una iteració
        Iteration++;
    }

    // Imprimim el nombre d'iteracions
    printf("Total d'iteracions: %i\n", Iteration - 1);
}


/*
* Funció main, crida als tres mètodes
*/
int main()
{

    // Creem les variables a usar
    int dimension = 1000000;
    double max_error = pow(10., -12.);
    double* vectorb;
    double* vector_solution;
    double* vector_solution_ant;
    clock_t start, end;
    double cpu_time_used;

    // Omplim el vector
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);

    // Mètode de Jacobi
    printf("\n\n********************\n");
    printf("*                  *\n");
    printf("* MÈTODE DE JACOBI *\n");
    printf("*                  *\n");
    printf("********************\n");

    // Iniciem el cronòmetre
    start = clock();

    // Cridem al mètode de Jacobi
    jacobi_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error);

    // Aturem el cronòmetre
    end = clock();

    // Calculem i imprimim el temps que s'ha trigat en fer totes les iteracions
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode de Jacobi: %f\n", cpu_time_used);


    // Tornem a inicialitzar els vectors 
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);

    // Mètode de Gauss-Seidel
    printf("\n\n**************************\n");
    printf("*                        *\n");
    printf("* MÈTODE DE GAUSS-SEIDEL *\n");
    printf("*                        *\n");
    printf("**************************\n");

    // Iniciem el cronòmetre
    start = clock();

    // Cridem al mètode de Gauss-Seidel
    gauss_seidel_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error);

    // Aturem el cronòmetre
    end = clock();

    // Calculem i imprimim el temps que s'ha trigat en fer totes les iteracions
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode de Gauss-Seidel: %f\n", cpu_time_used);


    // Tornem a omplir els vectors 
    vectorb = fill_vector(dimension);
    vector_solution = fill_vector_solution(dimension);
    vector_solution_ant = fill_vector_solution(dimension);

    // Mètode SOR
    printf("\n\n**************\n");
    printf("*            *\n");
    printf("* MÈTODE SOR *\n");
    printf("*            *\n");
    printf("**************\n");

    // Creem la millor omega (la manera com la hem trobat es troba al PDF i al codi comentat al final d'aquesta funció)
    double best_omega = 1.15;
    
    // Iniciem el cronòmetre
    start = clock();

    // Cridem al mètode SOR
    SOR_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error, best_omega);

    // Aturem el cronòmetre
    end = clock();

    // Imprimim l'omega usada
    printf("Omega usada: %.2f\n", best_omega);

    // Calculem i imprimim el temps que s'ha trigat en fer totes les iteracions
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode SOR: %f\n", cpu_time_used);


    // Codis per a trobar la millor omega


    // Codi per a provar omega entre 0 i 2 en escala 0.1
    /*
    double omega;
    start = clock();
    for(int k = 1; k < 20; k++){
        omega = k/10.;
        vectorb = fill_vector(dimension);
        vector_solution = fill_vector_solution(dimension);
        vector_solution_ant = fill_vector_solution(dimension);
        SOR_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error, omega);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode SOR: %f\n", cpu_time_used);
    */


    // Codi per a provar omega entre 1.1 i 1.2 en escala 0.01
    /*
    double omega;
    start = clock();
    for(int k = 0; k < 11; k++){
        omega = 1.1 + k/100.;
        vectorb = fill_vector(dimension);
        vector_solution = fill_vector_solution(dimension);
        vector_solution_ant = fill_vector_solution(dimension);
        SOR_method(vector_solution, vector_solution_ant, vectorb, dimension, max_error, omega);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Temps del mètode SOR: %f\n", cpu_time_used);
    */

    // Alliberem la memòria dels vectors
    free(vectorb);
    free(vector_solution);
    free(vector_solution_ant);

    // Tot ha anat bé, retornem un zero
    return 0;
}
