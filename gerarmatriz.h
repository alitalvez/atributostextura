#ifndef GEREMATRIZ_H
#define GEREMATRIZ_H

#include <iostream>

using namespace std;

class GerarMatriz
{
private:

public:
    ~GerarMatriz();
    unsigned short int* gerarMatriz(int tam, int nivel_cinza);
    unsigned short int* readFileMatriz(int tam);
    double* lerMatrizCoocorrencia(int tam);
};

#endif // GEREMATRIZ_H
