#ifndef ATCPU_H
#define ATCPU_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <stdio.h>
#include "tempo.h"

class atCpu {
private:
    int tamanho;
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

    double correlacao2(); // f3
    double somaEntropia2(); // f8
    double diferencaEntropia2(); // f11
    double mediaSoma2(); // f6
    double varianciaSoma2(); // f7
    double mediaDiferenca2();
    double varianciaDiferenca2(); // f10

public:

    atCpu(double *matriz, int tam) {
        this->matriz = matriz;
        this->tamanho = tam;
    }

    double max();
    double min();
    double media();
    double mediana();

    double variancia(); // f4
    double correlacao(); // f3
    double contraste(); // f2
    double entropia(); // f9
    double somaEntropia(); // f8
    double diferencaEntropia(); // f11
    double mdi(); // f5
    double mediaSoma(); // f6
    double varianciaSoma(); // f7
    double mediaDiferenca();
    double varianciaDiferenca(); // f10
    double energia(); // f1
    double MIC1();
    double MIC2();

    void CalcularAtributos();
};

#endif // ATCPU_H
