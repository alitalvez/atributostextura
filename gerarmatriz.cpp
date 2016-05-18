#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include "gerarmatriz.h"

GerarMatriz::~GerarMatriz() { }

unsigned short int* GerarMatriz::gerarMatriz(int tam, int nivel_cinza)
{
    srand(time(NULL));

    unsigned short int *matriz = new unsigned short int[tam*tam];

    if (nivel_cinza >= 2) {
        for(int i = 0; i < tam*tam; ++i) {
            matriz[i] = rand()%nivel_cinza;
        }
    } else {
        for(int i = 0; i < tam*tam; ++i)
            matriz[i] = 1;
    }

    return matriz;
}

unsigned short int* GerarMatriz::readFileMatriz(int tam)
{
    fstream myfile("in.txt", ios_base::in);

    unsigned short int *matriz = new unsigned short int[tam*tam];

    int i = 0;
    while (myfile >> matriz[i]) { i++; }

    myfile.close();

    return matriz;
}

double* GerarMatriz::lerMatrizCoocorrencia(int tam)
{
    fstream myfile("mc.txt", ios_base::in);

    double *matriz = new double[tam*tam];

    int i = 0;
    while (myfile >> matriz[i]) { i++; }

    myfile.close();

    return matriz;
}
