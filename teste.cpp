#include <iostream>
#include <cmath>

using namespace std;

double *matrizCoN = new double[16];

double px(int i)
{
    int tam  = 4;
    double px = 0.0;

//#pragma omp parallel for reduction(+:px)
    for(int k = 0; k < tam; ++k)
        px += matrizCoN[i * tam + k];
    return px;
}

double py(int j)
{
    int tam = 4;
    double py = 0.0;

//#pragma omp parallel for reduction(+:py)
    for(int k = 0; k < tam; ++k)
        py += matrizCoN[k * tam + j];
    return py;
}

double hxy()
{
    int tTotal = 16;
    double hxy = 0.0;
#pragma omp parallel for reduction(+:hxy)
    for(int k = 0; k < tTotal; ++k)
        if(matrizCoN[k])
            hxy += matrizCoN[k] * log(matrizCoN[k]);
    return hxy * -1;
}


double hxy2()
{
    int tam = 4;
    double hxy2 = 0.0;

    double rpx = 0.0;
    double rpy = 0.0;

#pragma omp parallel for reduction(+:hxy2)
    for(int k = 0; k < tam; ++k)
    {
        for(int h = 0; h < tam; ++h)
        {
            rpx = px(k);
            rpy = py(h);
            if(rpx && rpy)
            {
                hxy2 += rpx * rpy * log(rpx * rpy);
            }
        }
    }
    return hxy2 * -1;
}

double P_x_mais_y(const double * __restrict__ p, const int k, int tam)
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

double P_x_menos_y(const double * __restrict__ p, const int k, int tam)
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

double energia()
{
    int tTotal = 16;
    double energia = 0.0;

    //Tempo.Start();
#pragma omp parallel for simd reduction(+:energia)
    for(int i = 0; i < tTotal; ++i) {
        energia +=  matrizCoN[i] * matrizCoN[i];
    }

    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return energia;
}

double contraste()
{
    const int tam = 4;
    double contraste = 0.0;

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

double correlacao()
{
    double cor = 0.0;
    double var = 0.0;

    //Tempo.Start();

    double media = mediaH(matrizCoN, 4);

#pragma omp parallel for reduction(+:cor,var)
    for (int i = 0; i < 4; ++i) {
        double var1 = 0.0;
        double cor1 = 0.0;
        double tmp = (i - media) * (i - media);
#pragma omp simd reduction(+:var1,cor1)
        for (int j = 0; j < 4; ++j) {
            cor1 += (i - media) * (j - media) * matrizCoN[i * 4 + j];
            var1 += tmp * matrizCoN[i * 4 + j];
        }
        var += var1;
        cor += cor1;
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return cor/var;
}

double variancia()
{
    const int TAM = 4;
    double variancia = 0.0;

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

double mdi()
{
    const int tam = 4;
    double mdi = 0.0;

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

double mediaSoma()
{
    int tam = 4;
    int tTotal = (2*tam) - 2;
    double resultado = 0.0;

    //Tempo.Start();
#pragma omp parallel for reduction(+:resultado)
    for (int k = 0; k <= tTotal; ++k) {
        resultado += k * P_x_mais_y(matrizCoN, k, tam);
    }
    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return resultado;
}

double varianciaSoma()
{
    int tam = 4;
    int tTotal = (2*tam) - 2;
    double mediaS = 0.0;
    double varSoma = 0.0;

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

double somaEntropia()
{
    int tam = 4;
    const int tTotal = (2*tam) - 2;
    double somaEnt = 0.0;

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

double entropia()
{
    const int tTotal = 16;
    double entropia = 0.0;

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

double varianciaDiferenca()
{
    double * __restrict__ pxmy_temp = new double[4];
    double resultado = 0.0;

    //Tempo.Start();

/****/
#pragma omp parallel for
    for (int k = 0; k < 4; ++k) {
        pxmy_temp[k] = P_x_menos_y(matrizCoN, k, 4);
    }
/****/

#pragma omp parallel for simd reduction(+:resultado)
//#pragma omp simd reduction(+:resultado)
    for (int k=0; k < 4; ++k) {
        resultado += ( (k - k*pxmy_temp[k])*(k - k*pxmy_temp[k]) ) * pxmy_temp[k] ;
    }

    //Tempo.Stop();
    //Tempo.ElapsedSec();

    return resultado;
}

double diferencaEntropia()
{
    int tam = 4;
    double difEnt = 0.0;

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

double hx()
{
    int tam = 4;
    double hx = 0.0;
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

double hy()
{
    int tam = 4;
    double hy = 0.0;
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

double hxy1()
{
    int tam = 4;
    double hxy1 = 0.0;

#pragma omp parallel for reduction(+:hxy1)
    for(int k = 0; k < tam; ++k)
        for(int h = 0; h < tam; ++h)
            if(px(k) && py(h))
                hxy1 += matrizCoN[k * tam + h] * log(px(k) * py(h));
    return hxy1 * -1;
}

double medidasCorrelacao1()
{
    double rhxy = 0.0;
    double rhxy1 = 0.0;
    double mc = 0.0;

    rhxy = hxy();
    rhxy1 = hxy1();
    mc = (entropia() - rhxy1) / std::max(hx(), hy());

    return mc;
}

double medidasCorrelacao2()
{
    double mc = 0.0;
    double rhxy2 = 0.0;
    double rhxy = 0.0;

    rhxy2 = hxy2();
    rhxy = hxy();
    mc = std::sqrt(1 - std::exp(-2 * abs(rhxy2 - entropia())));

    return mc;
}


int main(int argc, char const *argv[]) {
    double result = 0.0;
    matrizCoN[0] = 0.25;
    matrizCoN[1] = 0;
    matrizCoN[2] = 0.083;
    matrizCoN[3] = 0;
    matrizCoN[4] = 0;
    matrizCoN[5] = 0.167;
    matrizCoN[6] = 0.083;
    matrizCoN[7] = 0;
    matrizCoN[8] = 0.083;
    matrizCoN[9] = 0.083;
    matrizCoN[10] = 0.083;
    matrizCoN[11] = 0.083;
    matrizCoN[12] = 0;
    matrizCoN[13] = 0;
    matrizCoN[14] = 0.083;
    matrizCoN[15] = 0;

    cout << "Energia: " << energia() << endl << endl;
    cout << "Contraste: " << contraste() << endl << endl;
    cout << "Entropia: " << entropia() << endl << endl;
    cout << "Variancia: " << variancia() << endl << endl;
    cout << "Correlacao: " << correlacao() << endl << endl;
    cout << "MDI: " << mdi() << endl << endl;
    cout << "Soma Media: " << mediaSoma() << endl << endl;
    cout << "Soma Variancia: " << varianciaSoma() << endl << endl;
    cout << "Soma Entropia: " << somaEntropia() << endl << endl;
    cout << "Dif Variancia: " << varianciaDiferenca() << endl << endl;
    cout << "Dif Entropia: " << diferencaEntropia() << endl << endl;
    cout << "Med Correlacao 1: " << medidasCorrelacao1() << endl << endl;
    cout << "Med Correlacao 2: " << medidasCorrelacao2() << endl << endl;

    return 0;
}
