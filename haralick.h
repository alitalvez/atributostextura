#ifndef HARALICK_H
#define HARALICK_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <stdio.h>

#include "globalvar.h"
#include "tempo.h"

class Haralick
{
private:
    int Ng;
    // Matriz de Coocorrencia
    int nLin, nCol;
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

    // Atributos

    double *matriz;

    CpuTimer Tempo;

    double *Pxpy;
    double *Pxmy;

    double Media;

    double Ener;
    double Con;
    double Cor;
    double var;
    double MDI;
    double MSoma;
    double SomaVar;
    double SomaEntr;
    double Entr;
    double DifVar;
    double DifEntr;

    void CalPxpy();
    void CalPxmy();

    double P_x_mais_y(const double * __restrict__ p, const int k, int tam);
    double P_x_menos_y(const double * __restrict__ p, const int k, int tam);

    double correlacao_(); // f3
    double somaEntropia_(); // f8
    double diferencaEntropia_(); // f11
    double mediaSoma_(); // f6
    double varianciaSoma_(); // f7
    double mediaDiferenca2();
    double varianciaDiferenca_(); // f10

public:

    Haralick(tImage* st_image) {
        this->matrizImg = (*st_image).vi_vector;
        this->nCol      = (*st_image).vi_coluna;
        this->nLin      = (*st_image).vi_linha;
        this->Ng        = pow(2,(*st_image).vi_bits);
    }

    ~Haralick() {

    }

    // Matriz de Coocorrencia

    void calcularMatrizCoN(double* matrizCoN, int distancia);

    // Atributos

    void atCpu(double *matriz, int tam) {
        this->matriz = matriz;
        this->Ng = tam;
    }

    double max();
    double min();
    double media();
    double mediana();

    double energia(); // f1
    double contraste(); // f2
    double correlacao(); // f3
    double variancia(); // f4
    double mdi(); // f5
    double mediaSoma(); // f6
    double varianciaSoma(); // f7
    double somaEntropia(); // f8
    double entropia(); // f9
    double varianciaDiferenca(); // f10
    double diferencaEntropia(); // f11
    double medidasCorrelacao1(); // f12
    double medidasCorrelacao2(); // f13

    double hx();
    double hy();
    double px(int);
    double py(int);
    double hxy();
    double hxy1();
    double hxy2();

    //double MIC1();
    //double MIC2();

    void CalcularAtributos();
};

#endif // HARALICK_H
