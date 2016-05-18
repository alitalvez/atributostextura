#ifndef MATRIZCOOCORRENCIA_H
#define MATRIZCOOCORRENCIA_H

//#include <iostream>
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <omp.h>

#include "globalvar.h"

using namespace std;

class MatrizCoocorrencia
{
private:
    int nLin, nCol;
    int Ng;
    int distancia;
    unsigned short *matrizImg;

    void mc0_(unsigned short* mIMG, int i, int j, int* matrizCoTmp);
    void mc45_(unsigned short* mIMG, int i, int j, int* matrizCoTmp);
    void mc90_(unsigned short* mIMG, int i, int j, int* matrizCoTmp);
    void mc135_(unsigned short* mIMG, int i, int j, int* matrizCoTmp);

    void mc0(int *matriz_freq);
    void mc45(int *matriz_freq);
    void mc90(int *matriz_freq);
    void mc135(int *matriz_freq);
    void normalizar(int*    __restrict__ matrizCo,
                    double* __restrict__ matrizCoN,
                    const int TAM_TOTAL);

public:
    MatrizCoocorrencia(tImage* st_image) {
        this->matrizImg = (*st_image).vi_vector;
        this->nCol      = (*st_image).vi_coluna;
        this->nLin      = (*st_image).vi_linha;
        this->Ng        = pow(2,(*st_image).vi_bits);
    };

    ~MatrizCoocorrencia() {
        //delete[] matriz;
    };

    void calcularMatrizCoN(double* matrizCoN, int distancia);
};

#endif // MATRIZCOOCORRENCIA_H
