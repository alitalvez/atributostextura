#include "matrizcoocorrencia.h"

void MatrizCoocorrencia::calcularMatrizCoN(double * __restrict__ matrizCoN, int distancia)
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
    int num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }
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

inline void MatrizCoocorrencia::mc0_(unsigned short * __restrict__ mIMG,
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

inline void MatrizCoocorrencia::mc45_(unsigned short * __restrict__ mIMG,
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
inline void MatrizCoocorrencia::mc90_(unsigned short * __restrict__ mIMG,
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

inline void MatrizCoocorrencia::mc135_(unsigned short * __restrict__ mIMG,
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


void MatrizCoocorrencia::mc0(int* matrizCoTmp)
{
	int tLin = this->nLin;
	int tCol = this->nCol;

    //int posicao = i + d;

#pragma omp parallel for
    for(int i = 0; i < tLin; i++)
        for(int j = 0; j < tCol; j++) {

            int k = j + distancia;

            if (k < nCol) {
                int z1 = matrizImg[ i*nCol + j];
                int z2 = matrizImg[ i*nCol + k];

                int ii = z1*Ng + z2;
//                int jj = z2*nivel_cinza + z1;

                #pragma omp atomic
                    matrizCoTmp[ii] += 1;
//                #pragma omp atomic
//                    matriz_freq[jj] += 1;
            }
        }
//}
    //return matriz_freq;
}

void MatrizCoocorrencia::mc45(int *matriz_freq)
{
	int tLin = this->nLin;
	int tCol = this->nCol;

    //int posicao = i - (n-1)*d;

#pragma omp parallel for
    for(int i = 0; i < tLin; i++)
        for(int j = 0; j < tCol; j++) {
            int ii = i - distancia;
            int jj = j + distancia;

            if (jj < nCol && ii >= 0) {
                int z1 = matrizImg[i*tCol + j];
                int z2 = matrizImg[ii*tCol + jj];

                ii = z1*Ng + z2;
//                jj = z2*nivel_cinza + z1;

                #pragma omp atomic
                    matriz_freq[ii] += 1;
//                #pragma omp atomic
//                    matriz_freq[jj] += 1;
            }
        }

    //return matriz_freq;
}

void MatrizCoocorrencia::mc90(int *matriz_freq)
{
	int tLin = this->nLin;
	int tCol = this->nCol;
	int tTotal = tLin * tCol;

    //int posicao = i - n*d;

#pragma omp parallel for
    for(int i = tCol * distancia; i < tTotal; ++i) {

        int posicao = i - tCol * distancia;

        if (posicao >= 0) {
            int z1 = matrizImg[i];
            int z2 = matrizImg[posicao];

            int ii = z1*Ng + z2;
            //int jj = z2*nivel_cinza + z1;

            #pragma omp atomic
                matriz_freq[ii] += 1;
            //#pragma omp atomic
                //matriz_freq[jj] += 1;
        }
    }

    //return matriz_freq;
}

void MatrizCoocorrencia::mc135(int *matriz_freq)
{

	int tLin = this->nLin;
	int tCol = this->nCol;
    //int posicao = i - (n+1)*d;

#pragma omp parallel for
    for(int i = 0; i < tLin; i++)
        for(int j = 0; j < tCol; j++) {
            int ii = i - distancia;
            int jj = j - distancia;

            if (jj >= 0 && ii >= 0) {

                int z1 = matrizImg[i*tCol + j];
                int z2 = matrizImg[ii*tCol + jj];

                ii = z1*Ng + z2;
                //jj = z2*nivel_cinza + z1;

                #pragma omp atomic
                    matriz_freq[ii] += 1;
                //#pragma omp atomic
                    //matriz_freq[jj] += 1;
            }
        }

    //return matriz_freq;
}

inline void MatrizCoocorrencia::normalizar(int*     __restrict__ matrizCo,
                                           double*  __restrict__ matrizCoN,
                                           const int TAM_TOTAL)
{
    int soma = 0;

//    #pragma omp parallel for simd reduction(+:soma)
    #pragma omp simd reduction(+:soma)
    for(int i = 0; i < TAM_TOTAL; ++i)
        soma += matrizCo[i];

//    printf("Soma CPU: %d \n", soma);

    double dsoma = (double)soma;
//    #pragma omp parallel for simd
    #pragma omp simd
    for(int i = 0; i < TAM_TOTAL; ++i)
        matrizCoN[i] = 1.0*matrizCo[i] / dsoma;

}
