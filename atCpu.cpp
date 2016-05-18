#include "atCpu.h"

/* Protótipo */
//static int compareMyType (const void * a, const void * b);

//static double Px(const double *p, const int i, int tam);
//static double Py(const double *p, const int j, int tam);
//static double P_x_mais_y(const double *p, const int k, int tam);
//static double P_x_menos_y(const double *p, const int k, int tam);
static double mediaH(const double * __restrict__ p, int tam);

/* Definições */

/*
int compareMyType (const void * a, const void * b)
{
    if ( *(double*)a <  *(double*)b ) return -1;
    if ( *(double*)a == *(double*)b ) return 0;

    return 1;
}

double Px(const double *p, const int i, int tam)
{
    double result = 0.0;

    for (int j=0; j < tam; ++j) {
        result += p[i * tam + j];
    }
    return result;
}

double Py(const double *p, const int j, int tam)
{
    double result = 0.0;

    for (int i=0; i < tam; ++i) {
        result += p[i * tam + j];
    }
    return result;
}
*/

/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

void atCpu::CalPxpy()
{
    const double * __restrict__ p = this->matriz;
    const int tam = this->tamanho;
    int tTotal = (2*tam) - 2;

    Pxpy = new double[tTotal];

//#pragma omp parallel for
    for (int k = 0; k < tTotal; ++k) {
        double soma = 0.0;
        for (int i=0; i < tam; ++i) {
            double s = 0.0;
#pragma omp simd reduction(+:s)
            for (int j=0; j < tam; ++j) {
                if ( (i + j) == k)
                    s += p[i * tam + j];
            }
            soma += s;
        }
        Pxpy[k] = soma;
    }

}

void atCpu::CalPxmy()
{
    const double * __restrict__ p = this->matriz;
    const int tam = this->tamanho;
    int tTotal = tam;

    Pxmy = new double[tTotal];

//#pragma omp parallel for
    for (int k = 0; k < tTotal; ++k) {
        double soma = 0.0;
        for (int i=0; i < tam; ++i) {
            double s = 0.0;
#pragma omp simd reduction(+:s)
            for (int j=0; j < tam; ++j) {
                if ( abs(i - j) == k)
                    s += p[i * tam + j];
            }
            soma += s;
        }
        Pxmy[k] = soma;
    }

}

void atCpu::CalcularAtributos()
{

//    CpuTimer T;

//    T.Start();

    atCpu::CalPxmy();
    atCpu::CalPxpy();
    Media = mediaH(this->matriz, this->tamanho);

    Ener = atCpu::energia(); //f1
//    printf("Ener:     %lf\n", Ener);
    MDI  = atCpu::mdi(); //f5
//    printf("MDI:      %lf\n", MDI);
    Entr = atCpu::entropia(); //f9
//    printf("Entr:     %lf\n", Entr);
    var  = atCpu::variancia(); // f4
//    printf("Var:      %lf\n", var);
    Con  = atCpu::contraste(); // f2
//    printf("Con:      %lf\n", Con);
    Cor  = atCpu::correlacao2(); // f3
//    printf("Cor:      %lf\n", Cor);
    MSoma = atCpu::mediaSoma2(); // f6
//    printf("MSoma:    %lf\n", MSoma);
    SomaVar  = atCpu::varianciaSoma2(); // f7
//    printf("SomaVar:  %lf\n", SomaVar);
    SomaEntr = atCpu::somaEntropia2(); // f8
//    printf("SomaEntr: %lf\n", SomaEntr);
    DifVar   = atCpu::varianciaDiferenca2(); // f10
//    printf("DifVar:   %lf\n", DifVar);
    DifEntr  = atCpu::diferencaEntropia2(); // f11
//    printf("DifEntr:  %lf\n", DifEntr);

//    T.Stop();
//    T.ElapsedSec();

}

/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

double atCpu::P_x_mais_y(const double * __restrict__ p, const int k, int tam)
{
    double result = 0.0;

    for (int i=0; i < tam; ++i) {
        #pragma omp simd reduction(+:result)
        for (int j=0; j < tam; ++j)
            if ( (i + j) == k)
                result += p[i * tam + j];
    }

    return result;
}

double atCpu::P_x_menos_y(const double * __restrict__ p, const int k, int tam)
{
    double result = 0.0;

    for (int i=0; i < tam; ++i) {
    #pragma omp simd reduction(+:result)
        for (int j=0; j < tam; ++j)
            if ( std::abs(i - j) == k)
                result += p[i * tam + j];
    }

    return result;
}


double mediaH(const double * __restrict__ p, int tam)
{
    double soma = 0.0;

//#pragma omp parallel for reduction(+:soma)
    for (int i = 0; i < tam; ++i) { 
        double soma1 = 0.0;
#pragma omp simd reduction(+:soma1)
        for (int j = 0; j < tam; ++j)
            soma1 += i * p[i * tam + j];
        soma += soma1;
    }

    return soma;
}

/*
 * Atributos
 */

/*
double atCpu::max()
{
    int TAM_M = tamanho * tamanho;
    double valor = 0.0;

    for(int i=0; i < TAM_M; i++) {
        valor = fmax(matriz[i], valor);
    }

    return valor;
}

double atCpu::min()
{
    int TAM_M = tamanho * tamanho;
    double valor = 1.0;

    for(int i=0; i < TAM_M; i++) {
        valor = fmin(matriz[i], valor);
    }

    return valor;
}

double atCpu::media()
{
    int TAM_M = tamanho * tamanho;
    double media = 0.0;

    for (int i = 0; i < TAM_M; i++) {
        media += matriz[i];
    }

    //std::cout << "tm: " << media << std::endl;

    media = media / (TAM_M);

    return media;
}

double atCpu::mediana()
{
    int tTotal = tamanho * tamanho;
    double mediana;
    double *array =  new double[tTotal];

#pragma omp parallel for //simd
    for (int i=0; i < tTotal; i++) {
        array[i] = matriz[i];
    }

    std::qsort(array, (size_t) tTotal, sizeof(double), compareMyType );

    unsigned pos = (tTotal)/2;
    mediana = (array[ pos - 1 ] + array[ pos ]) / 2;

    delete [] array;

    return mediana;
}
*/

/*
 * Haralick
 */
double atCpu::variancia()
{
    const int tam = this->tamanho;
    double variancia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//    Tempo.Start();
//    double media = this->Media; //mediaH(matrizCoN, tam);
    double media = mediaH(matrizCoN, tam);

//#pragma omp parallel for reduction(+:variancia)
    for (int i = 0; i < tam; i++) {
        double var = 0.0;
        double tmp = (i - media) * (i - media);
#pragma omp simd reduction(+:var)
        for(int j = 0; j < tam; j++)
            var += tmp * matrizCoN[i * tam + j];
        variancia += var;
    }
//    Tempo.Stop();
//    Tempo.ElapsedSec();

    return variancia;
}

double atCpu::correlacao2()
{
    const int tam = this->tamanho;
    double cor = 0.0;
    double var = this->var; //0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//	Tempo.Start();

    double media = this->Media;

//#pragma omp parallel for reduction(+:cor)
    for (int i = 0; i < tam; ++i) {
        double cor1 = 0.0;
#pragma omp simd reduction(+:cor1)
        for (int j = 0; j < tam; ++j) {
            cor1 += (i*j) * (i - media) * (j - media) * matrizCoN[i * tam + j];
        }
        cor += cor1;
    }
//	Tempo.Stop();
//    Tempo.ElapsedSec();

    return cor/var;
}

double atCpu::correlacao()
{
    const int tam = this->tamanho;
    double cor = 0.0;
    double var = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();
    double media = mediaH(matrizCoN, tam);
//	Tempo.Stop();
//    Tempo.ElapsedSec();

#pragma omp parallel for reduction(+:cor,var)
    for (int i = 0; i < tam; ++i) {
        double var1 = 0.0;
        double cor1 = 0.0;
        double tmp = (i - media) * (i - media);
#pragma omp simd reduction(+:var1,cor1)
        for (int j = 0; j < tam; ++j) {
            cor1 += (i*j) * (i - media) * (j - media) * matrizCoN[i * tam + j];
            var1 += tmp * matrizCoN[i * tam + j];
        }
        var += var1;
        cor += cor1;
	}
    Tempo.Stop();
    Tempo.ElapsedSec();
	
    return cor/var;
}

double atCpu::contraste()
{
    const int tam = this->tamanho;
    double contraste = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//    Tempo.Start();
//#pragma omp parallel for reduction(+:contraste)
    for (int i = 0; i < tam; ++i) {
        double con = 0.0;
#pragma omp simd reduction(+:con)
        for(int j = 0; j < tam; ++j) {
            con += ( std::abs(i - j) * std::abs(i - j)) * matrizCoN[i * tam + j];
        }
        contraste += con;
    }
//    Tempo.Stop();
//    Tempo.ElapsedSec();

    return contraste;
}

double atCpu::entropia()
{
    const int tTotal = tamanho * tamanho;
    double entropia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//    Tempo.Start();
//#pragma omp parallel for simd reduction(+:entropia)
#pragma omp simd reduction(+:entropia)
    for (int i = 0; i < tTotal; ++i) {
        double valor = matrizCoN[i];
        entropia += (valor > 0.0) ? valor * log2(valor) : 0.0;
    }
//    Tempo.Stop();
//    Tempo.ElapsedSec();

    return entropia * (-1.0);
}


double atCpu::somaEntropia2()
{
    int tam = this->tamanho;
    const int tTotal = (2*tam) - 2;
    double somaEnt = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

//#pragma omp parallel for simd reduction(+:somaEnt)
#pragma omp simd reduction(+:somaEnt)
    for (int k = 0; k <= tTotal; ++k) {
        if (pxpy[k] > 0.0)
            somaEnt += pxpy[k] * log2(pxpy[k]);
    }

    return (-1.0) * somaEnt;
}

double atCpu::somaEntropia()
{
    int tam = this->tamanho;
    const int tTotal = (2*tam) - 2;
    double somaEnt = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();
#pragma omp parallel for reduction(+:somaEnt)
    for (int k = 0; k <= tTotal; ++k) {
        double tmp = P_x_mais_y(matrizCoN, k, tam);
        if (tmp > 0.0)
            somaEnt += tmp * log2(tmp);
    }
    Tempo.Stop();
    Tempo.ElapsedSec();

    return (-1.0) * somaEnt;
}


double atCpu::diferencaEntropia2()
{
    int tam = this->tamanho;
    double difEnt = 0.0;
    double * __restrict__ pxmy = this->Pxmy;

//#pragma omp parallel for simd reduction(+:difEnt)
#pragma omp simd reduction(+:difEnt)
    for (int k = 0; k < tam; ++k) {
        if (pxmy[k] > 0.0)
            difEnt += pxmy[k] * log2(pxmy[k]);
    }

    return (-1.0) * difEnt;
}

double atCpu::diferencaEntropia()
{
    int tam = this->tamanho;
    double difEnt = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

	Tempo.Start(); 

#pragma omp parallel for reduction(+:difEnt)
    for (int k = 0; k < tam; ++k) {
        double tmp = P_x_menos_y(matrizCoN, k, tam);
        if (tmp > 0.0)
            difEnt += tmp * log2(tmp);
    }
    Tempo.Stop();
    Tempo.ElapsedSec();

    return (-1.0) * difEnt;
}

double atCpu::mdi()
{
    const int tam = this->tamanho;
    double mdi = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//    Tempo.Start();
//#pragma omp parallel for reduction(+:mdi)
    for(int i = 0; i < tam; ++i) {
        double mdi1 = 0.0;
#pragma omp simd reduction(+:mdi1)
        for(int j = 0; j < tam; ++j) {
            mdi1 += ( matrizCoN[i * tam + j] / (1 + ((i - j)*(i - j)) ) );
        }
        mdi += mdi1;
    }
//    Tempo.Stop();
//    Tempo.ElapsedSec();

    return mdi;
}

double atCpu::mediaSoma2()
{
    int tam = this->tamanho;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

//#pragma omp parallel for simd reduction(+:resultado)
#pragma omp simd reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * pxpy[k];
    }

    return resultado;
}

double atCpu::mediaSoma()
{
    int tam = this->tamanho;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();
#pragma omp parallel for reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * P_x_mais_y(matrizCoN, k, tam);
    }
    Tempo.Stop();
    Tempo.ElapsedSec();

    return resultado;
}

double atCpu::varianciaSoma2()
{
    int tam = this->tamanho;
    int tTotal = (2*tam) - 2;
    double varSoma = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

    double mediaS = MSoma;

//#pragma omp parallel for simd reduction(+:varSoma)
#pragma omp simd reduction(+:varSoma)
    for (int k = 0; k <= tTotal; ++k) {
        varSoma += (k - mediaS) * (k - mediaS) * pxpy[k];
    }

    return varSoma;
}

double atCpu::varianciaSoma()
{
    int tam = this->tamanho;
    int tTotal = (2*tam) - 2;
    double mediaS = 0.0;
    double varSoma = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();

    mediaS = mediaSoma();

#pragma omp parallel for reduction(+:varSoma)
    for (int k = 0; k <= tTotal; ++k) {
        varSoma += (k - mediaS) * (k - mediaS) * P_x_mais_y(matrizCoN, k, tam);
    }
    Tempo.Stop();
    Tempo.ElapsedSec();

    return varSoma;
}

/*
double atCpu::mediaDiferenca()
{
    int tam = this->tamanho;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

#pragma omp parallel for simd reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * P_x_menos_y(matrizCoN, k, tam);
    }

    return resultado;
}
*/

double atCpu::varianciaDiferenca2()
{
    double resultado = 0.0;
    const double * __restrict__ pxmy = this->Pxmy;

#pragma omp parallel for simd reduction(+:resultado)
    for (int k=0; k < tamanho; ++k) {
        resultado += ( (k - k*pxmy[k]) * (k - k*pxmy[k]) ) * pxmy[k];
    }

    return resultado;
}

double atCpu::varianciaDiferenca()
{
    double * __restrict__ pxmy_temp = new double[tamanho];
    double resultado = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();

/****/
#pragma omp parallel for
    for (int k = 0; k < tamanho; ++k) {
        pxmy_temp[k] = P_x_menos_y(matrizCoN, k, tamanho);
    }
/****/

#pragma omp parallel for simd reduction(+:resultado)
//#pragma omp simd reduction(+:resultado)
    for (int k=0; k < tamanho; ++k) {
        resultado += ( (k - k*pxmy_temp[k])*(k - k*pxmy_temp[k]) ) * pxmy_temp[k] ;
    }

    Tempo.Stop();
    Tempo.ElapsedSec();

    return resultado;
}

double atCpu::energia()
{
    int tTotal = this->tamanho * this->tamanho;
    double energia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

//    Tempo.Start();
#pragma omp parallel for simd reduction(+:energia)
    for(int i = 0; i < tTotal; ++i) {
        energia +=  matrizCoN[i] * matrizCoN[i];
    }
//    Tempo.Stop();
//    Tempo.ElapsedSec();

    return energia;
}
// */
