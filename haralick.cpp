#include "haralick.h"

void Haralick::calcularMatrizCoN(double * __restrict__ matrizCoN, int distancia)
{
    int imgLin = this->nLin;
    int imgCol = this->nCol;
    unsigned short * __restrict__ mIMG = this->matrizImg;

    const int N = this->Ng;

//    int * __restrict__ matrizCoTmp  = new int[N*N];
    int * __restrict__ matrizCoF = new int[N*N];

    this->distancia = distancia;

//    memset(matrizCoTmp, 0, N*N);

///////////////////////////////////
    int num_threads = 0;
    #pragma omp parallel
        num_threads = omp_get_num_threads();
    printf("Threads: %d\n", num_threads);

    int **mc;
    mc = new int*[num_threads];
    for(int i = 0; i < num_threads; i++) {
        mc[i] = new int[N*N];
    }

/////////////////////////////////////
    #pragma omp parallel
    {
        //int tid = 0;
        int tid = omp_get_thread_num();
        #pragma omp for //schedule(dynamic,1)
        for(int i = 0; i < imgLin; ++i) {
            for(int j = 0; j < imgCol; ++j) {

                mc0_    (mIMG, i, j, mc[tid]);
                mc45_   (mIMG, i, j, mc[tid]);
                mc90_   (mIMG, i, j, mc[tid]);
                mc135_  (mIMG, i, j, mc[tid]);

    //            mc0_(mIMG, i, j, matrizCoTmp); // 7687830
    //            mc45_(mIMG, i, j, matrizCoTmp); // 7686000
    //            mc90_(mIMG, i, j, matrizCoTmp); // 7690200
    //            mc135_(mIMG, i, j, matrizCoTmp); // 7686000
            }
        }
    }
    ///*
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            for(int k = 1; k < num_threads; k++)
                mc[0][i*N + j] += mc[k][i*N + j];
                //matrizCoTmp[i*N + j] += mc[k][i*N + j];
    }
    //*/

    // Transposta
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        #pragma omp simd
        for (int j = 0; j < N; ++j)
            matrizCoF[i*N + j] = mc[0][i*N + j] + mc[0][j*N + i];
            //matrizCoF[i*N + j] = matrizCoTmp[i*N + j] + matrizCoTmp[j*N + i];
    }

    // Normalizacao

    const int TAM_TOTAL = N*N;
    int soma = 0;

    //#pragma omp simd reduction(+:soma)
    #pragma omp parallel for simd reduction(+:soma)
    for(int i = 0; i < TAM_TOTAL; ++i)
        soma += matrizCoF[i];

    printf("Soma CPU: %d \n", soma);

    //#pragma omp simd
    #pragma omp parallel for simd
    for(int i = 0; i < TAM_TOTAL; ++i)
        matrizCoN[i] = 1.0*matrizCoF[i] / soma;
    //normalizar(matrizCoF, matrizCoN, N*N);

//    delete [] matrizCoTmp;
    delete [] matrizCoF;
}

inline void Haralick::mc0_(unsigned short * __restrict__ mIMG,
                                     int i, int j,
                                     int* matrizCoTmp)
{
    int k = j + distancia;

    if (k < nCol) {
        int z1 = mIMG[ i*nCol + j];
        int z2 = mIMG[ i*nCol + k];

        int pos = z1*Ng + z2;

//        #pragma omp atomic
            matrizCoTmp[pos] += 1;
    }
}

inline void Haralick::mc45_(unsigned short * __restrict__ mIMG,
                                     int i, int j,
                                     int* matrizCoTmp)
{
    int ii = i - distancia;
    int jj = j + distancia;

    if (jj < nCol && ii >= 0) {
        int z1 = mIMG[i*nCol + j];
        int z2 = mIMG[ii*nCol + jj];

        int pos = z1*Ng + z2;

//        #pragma omp atomic
            matrizCoTmp[pos] += 1;
    }

}
inline void Haralick::mc90_(unsigned short * __restrict__ mIMG,
                                     int i, int j,
                                     int* matrizCoTmp)
{
    int ii = i - distancia;

    if (ii >= 0 && ii < nLin ) {
        int z1 = mIMG[i*nCol + j];
        int z2 = mIMG[ii*nCol + j];

        int pos = z1*Ng + z2;

//        #pragma omp atomic
            matrizCoTmp[pos] += 1;
    }
}

inline void Haralick::mc135_(unsigned short * __restrict__ mIMG,
                                     int i, int j,
                                     int* matrizCoTmp)
{
    int ii = i - distancia;
    int jj = j - distancia;

    if (jj >= 0 && ii >= 0) {

        int z1 = mIMG[i*nCol + j];
        int z2 = mIMG[ii*nCol + jj];

        int pos = z1*Ng + z2;

//        #pragma omp atomic
            matrizCoTmp[pos] += 1;
    }
}


/*
 * ****************************************************************************
 * ****************************************************************************
 * ****************************************************************************
 */


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

void Haralick::CalPxpy()
{
    const double * __restrict__ p = this->matriz;
    const int tam = this->Ng;
    int tTotal = (2*tam) - 2;

    Pxpy = new double[tTotal];

#pragma omp parallel for
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

void Haralick::CalPxmy()
{
    const double * __restrict__ p = this->matriz;
    const int tam = this->Ng;
    int tTotal = tam;

    Pxmy = new double[tTotal];

#pragma omp parallel for
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

void Haralick::CalcularAtributos()
{

//    CpuTimer T;
//    T.Start();

    Haralick::CalPxmy();
    Haralick::CalPxpy();
    Media = mediaH(this->matriz, this->Ng);

    Ener = Haralick::energia(); //f1
//    printf("Ener:     %lf\n", Ener);
    Con  = Haralick::contraste(); // f2
//    printf("Con:      %lf\n", Con);
    Cor  = Haralick::correlacao_(); // f3
//    printf("Cor:      %lf\n", Cor);
    var  = Haralick::variancia(); // f4
//    printf("Var:      %lf\n", var);
    MDI  = Haralick::mdi(); //f5
//    printf("MDI:      %lf\n", MDI);
    MSoma = Haralick::mediaSoma_(); // f6
//    printf("MSoma:    %lf\n", MSoma);
    SomaVar  = Haralick::varianciaSoma_(); // f7
//    printf("SomaVar:  %lf\n", SomaVar);
    SomaEntr = Haralick::somaEntropia_(); // f8
//    printf("SomaEntr: %lf\n", SomaEntr);
    Entr = Haralick::entropia(); //f9
//    printf("Entr:     %lf\n", Entr);
    DifVar   = Haralick::varianciaDiferenca_(); // f10
//    printf("DifVar:   %lf\n", DifVar);
    DifEntr  = Haralick::diferencaEntropia_(); // f11
//    printf("DifEntr:  %lf\n", DifEntr);

//    T.Stop();
//    T.ElapsedSec();

}

/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

double Haralick::P_x_mais_y(const double * __restrict__ p, const int k, int tam)
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

double Haralick::P_x_menos_y(const double * __restrict__ p, const int k, int tam)
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

#pragma omp parallel for reduction(+:soma)
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
double Haralick::max()
{
    int TAM_M = tamanho * tamanho;
    double valor = 0.0;

    for(int i=0; i < TAM_M; i++) {
        valor = fmax(matriz[i], valor);
    }

    return valor;
}

double Haralick::min()
{
    int TAM_M = tamanho * tamanho;
    double valor = 1.0;

    for(int i=0; i < TAM_M; i++) {
        valor = fmin(matriz[i], valor);
    }

    return valor;
}

double Haralick::media()
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

double Haralick::mediana()
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

double Haralick::energia()
{
    int tTotal = this->Ng * this->Ng;
    double energia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for simd reduction(+:energia)
    for(int i = 0; i < tTotal; ++i) {
        energia +=  matrizCoN[i] * matrizCoN[i];
    }

    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return energia;
}

double Haralick::contraste()
{
    const int tam = this->Ng;
    double contraste = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for reduction(+:contraste)
    for (int i = 0; i < tam; ++i) {
        double con = 0.0;
#pragma omp simd reduction(+:con)
        for(int j = 0; j < tam; ++j) {
            con += ( std::abs(i - j) * std::abs(i - j)) * matrizCoN[i * tam + j];
        }
        contraste += con;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return contraste;
}

double Haralick::correlacao_()
{
    const int tam = this->Ng;
    double cor = 0.0;
    double var = this->var; //0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    Tempo.Start();

    double media = this->Media;

#pragma omp parallel for reduction(+:cor)
    for (int i = 0; i < tam; ++i) {
        double cor1 = 0.0;
#pragma omp simd reduction(+:cor1)
        for (int j = 0; j < tam; ++j) {
            cor1 += (i*j) * (i - media) * (j - media) * matrizCoN[i * tam + j];
        }
        cor += cor1;
    }
    Tempo.Stop();
    Tempo.ElapsedSec();

    return cor/var;
}

double Haralick::correlacao()
{
    const int tam = this->Ng;
    double cor = 0.0;
    double var = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();

    double media = mediaH(matrizCoN, tam);

#pragma omp parallel for reduction(+:cor,var)
    for (int i = 0; i < tam; ++i) {
        double var1 = 0.0;
        double cor1 = 0.0;
        double tmp = (i - media) * (i - media);
#pragma omp simd reduction(+:var1,cor1)
        for (int j = 0; j < tam; ++j) {
            cor1 += (i - media) * (j - media) * matrizCoN[i * tam + j];
            var1 += tmp * matrizCoN[i * tam + j];
        }
        var += var1;
        cor += cor1;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return cor/var;
}


double Haralick::variancia()
{
    const int TAM = this->Ng;
    double variancia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();

    //double media = this->Media; //mediaH(matrizCoN, tam);
    double media = mediaH(matrizCoN, TAM);

#pragma omp parallel for reduction(+:variancia)
    for (int i = 0; i < TAM; i++) {
        double var = 0.0;
        double tmp = (i - media) * (i - media);
#pragma omp simd reduction(+:var)
        for(int j = 0; j < TAM; j++)
            var += tmp * matrizCoN[i * TAM + j];
        variancia += var;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return variancia;
}

double Haralick::mdi()
{
    const int tam = this->Ng;
    double mdi = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for reduction(+:mdi)
    for(int i = 0; i < tam; ++i) {
        double mdi1 = 0.0;
#pragma omp simd reduction(+:mdi1)
        for(int j = 0; j < tam; ++j) {
            mdi1 += ( matrizCoN[i * tam + j] / (1 + ((i - j)*(i - j)) ) );
        }
        mdi += mdi1;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return mdi;
}

double Haralick::mediaSoma_()
{
    int tam = this->Ng;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

#pragma omp parallel for simd reduction(+:resultado)
//#pragma omp simd reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * pxpy[k];
    }

    return resultado;
}

double Haralick::mediaSoma()
{
    int tam = this->Ng;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * P_x_mais_y(matrizCoN, k, tam);
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return resultado;
}

double Haralick::varianciaSoma_()
{
    int tam = this->Ng;
    int tTotal = (2*tam) - 2;
    double varSoma = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

    double mediaS = MSoma;

#pragma omp parallel for simd reduction(+:varSoma)
//#pragma omp simd reduction(+:varSoma)
    for (int k = 0; k <= tTotal; ++k) {
        varSoma += (k - mediaS) * (k - mediaS) * pxpy[k];
    }

    return varSoma;
}

double Haralick::varianciaSoma()
{
    int tam = this->Ng;
    int tTotal = (2*tam) - 2;
    double mediaS = 0.0;
    double varSoma = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();

    mediaS = mediaSoma();

#pragma omp parallel for reduction(+:varSoma)
    for (int k = 0; k <= tTotal; ++k) {
        varSoma += (k - mediaS) * (k - mediaS) * P_x_mais_y(matrizCoN, k, tam);
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return varSoma;
}


double Haralick::somaEntropia_()
{
    int tam = this->Ng;
    const int tTotal = (2*tam) - 2;
    double somaEnt = 0.0;
    double * __restrict__ pxpy = this->Pxpy;

#pragma omp parallel for simd reduction(+:somaEnt)
//#pragma omp simd reduction(+:somaEnt)
    for (int k = 0; k <= tTotal; ++k) {
        if (pxpy[k] > 0.0)
            somaEnt += pxpy[k] * log(pxpy[k]);
    }

    return (-1.0) * somaEnt;
}

double Haralick::somaEntropia()
{
    int tam = this->Ng;
    const int tTotal = (2*tam) - 2;
    double somaEnt = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for reduction(+:somaEnt)
    for (int k = 0; k <= tTotal; ++k) {
        double tmp = P_x_mais_y(matrizCoN, k, tam);
        if (tmp > 0.0)
            somaEnt += tmp * log(tmp);
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return (-1.0) * somaEnt;
}

double Haralick::entropia()
{
    const int tTotal = Ng * Ng;
    double entropia = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();
#pragma omp parallel for simd reduction(+:entropia)
//#pragma omp simd reduction(+:entropia)
    for (int i = 0; i < tTotal; ++i) {
        double valor = matrizCoN[i];
        entropia += (valor > 0.0) ? valor * log(valor) : 0.0;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return entropia * (-1.0);
}

double Haralick::varianciaDiferenca_()
{
    double resultado = 0.0;
    const double * __restrict__ pxmy = this->Pxmy;

#pragma omp parallel for simd reduction(+:resultado)
    for (int k=0; k < Ng; ++k) {
        resultado += ( (k - k*pxmy[k]) * (k - k*pxmy[k]) ) * pxmy[k];
    }

    return resultado;
}

double Haralick::varianciaDiferenca()
{
    double * __restrict__ pxmy_temp = new double[Ng];
    double resultado = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();

/****/
#pragma omp parallel for
    for (int k = 0; k < Ng; ++k) {
        pxmy_temp[k] = P_x_menos_y(matrizCoN, k, Ng);
    }
/****/

#pragma omp parallel for simd reduction(+:resultado)
//#pragma omp simd reduction(+:resultado)
    for (int k=0; k < Ng; ++k) {
        resultado += ( (k - k*pxmy_temp[k])*(k - k*pxmy_temp[k]) ) * pxmy_temp[k] ;
    }

    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return resultado;
}

double Haralick::diferencaEntropia_()
{
    int tam = this->Ng;
    double difEnt = 0.0;
    double * __restrict__ pxmy = this->Pxmy;

#pragma omp parallel for simd reduction(+:difEnt)
//#pragma omp simd reduction(+:difEnt)
    for (int k = 0; k < tam; ++k) {
        if (pxmy[k] > 0.0)
            difEnt += pxmy[k] * log(pxmy[k]);
    }

    return (-1.0) * difEnt;
}

double Haralick::diferencaEntropia()
{
    int tam = this->Ng;
    double difEnt = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

    //Tempo.Start();

#pragma omp parallel for reduction(+:difEnt)
    for (int k = 0; k < tam; ++k) {
        double tmp = P_x_menos_y(matrizCoN, k, tam);
        if (tmp > 0.0)
            difEnt += tmp * log(tmp);
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return (-1.0) * difEnt;
}

double Haralick::px(int i)
{
    int tam  = this->Ng;
    double px = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

#pragma omp parallel for reduction(+:px)
    for(int k = 0; k < tam; ++k)
        px += matrizCoN[i * tam + k];
    return px;
}

double Haralick::py(int j)
{
    int tam = this->Ng;
    double py = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

#pragma omp parallel for reduction(+:py)
    for(int k = 0; k < tam; ++k)
        py += matrizCoN[k * tam + j];
    return py;
}

double Haralick::hx()
{
    int tam = this->Ng;
    double hx = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;
    double rpx = 0.0;

#pragma omp parallel for reduction(+:hx) private(rpx)
    for(int h = 0; h < tam; ++h)
    {
        rpx = px(h);
        if(rpx)
            hx += rpx * log(rpx);
    }

    return hx * -1;
}

double Haralick::hy()
{
    int tam = this->Ng;
    double hy = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;
    double rpy = 0.0;

#pragma omp parallel for reduction(+:hy) private(rpy)
    for(int h = 0; h < tam; ++h)
    {
        rpy = py(h);
        if(rpy)
            hy += rpy * log(rpy);
    }
    return hy * -1;
}

double Haralick::hxy()
{
    int tTotal = this->Ng * this->Ng;
    double hxy = 0.0;
    double rlog = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;
#pragma omp parallel for reduction(+:hxy) private(rlog)
    for(int k = 0; k < tTotal; ++k)
        if(matrizCoN[k])
        {
            rlog = log(matrizCoN[k]);
            hxy += matrizCoN[k] * rlog;
        }
    return hxy * -1;
}

double Haralick::hxy1()
{
    int tam = this->Ng;
    double hxy1 = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;

#pragma omp parallel for reduction(+:hxy1)
    for(int k = 0; k < tam; ++k)
        for(int h = 0; h < tam; ++h)
            if(px(k) && py(h))
                hxy1 += matrizCoN[k * tam + h] * log(px(k) * py(h));
    return hxy1 * -1;
}

double Haralick::hxy2()
{
    int tam = this->Ng;
    double hxy2 = 0.0;
    const double * __restrict__ matrizCoN = this->matriz;
    double rpx = 0.0;
    double rpy = 0.0;
    double rlog = 0.0;

#pragma omp parallel for reduction(+:hxy2) private(rpx, rpy, rlog)
    for(int k = 0; k < tam; ++k)
    {
        for(int h = 0; h < tam; ++h)
        {
            rpx = px(k);
            rpy = py(h);
            if(rpx && rpy)
            {
                rlog = log(rpx * rpy);
                hxy2 += rpx * rpy * rlog;
            }
        }
    }
    return hxy2 * -1;
}

double Haralick::medidasCorrelacao1()
{
    double rhxy = 0.0;
    double rhxy1 = 0.0;
    double mc = 0.0;

    rhxy = hxy();
    rhxy1 = hxy1();
    mc = (entropia() - rhxy1) / std::max(hx(), hy());

    return mc;
}

double Haralick::medidasCorrelacao2()
{
    double mc = 0.0;
    double rhxy2 = 0.0;
    double rhxy = 0.0;

    rhxy2 = hxy2();
    rhxy = hxy();

    mc = std::sqrt(1 - std::exp(-2 * abs(rhxy2 - entropia())));

    return mc;
}
