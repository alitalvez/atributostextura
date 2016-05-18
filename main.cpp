//g++ -fopenmp -O3 -mavx libHara.o libRead.o libGera.o libPy.o libAt.o main.cpp -o tentativa.o

//#include <QCoreApplication>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "readimage.h"
//#include "gerarmatriz.h"
//#include "matrizcoocorrencia.h"
#include "haralick.h"
#include "tempo.h"
#include "utilidades.h"

int main(int argc, char *argv[])
{
    char filename[32];
    int coluna;
    int linha;
    int Ng = 4096;
    int distancia = 1;
    int op = 0;

    double startTime;
    double stopTime;

    //QCoreApplication a(argc, argv);

    tImage st_image;
    /*
    if (argc == 1) {
        coluna = 3000;
        linha = 3000;

        GerarMatriz g_imagem;
        st_image.vi_coluna = coluna;
        st_image.vi_linha = linha;
        //st_image.vi_vector = g_imagem.gerarMatriz(width, nivelCinza);
        st_image.vi_vector = g_imagem.readFileMatriz(coluna);
        st_image.vi_bits = 12;

    } else
    */

    if (argc == 7) {
        strcpy(filename, argv[1]);
        linha       = atoi(argv[2]);
        coluna      = atoi(argv[3]);
        Ng          = (int) pow(2, atoi(argv[4])); // Nivel de Cinza
        distancia   = atoi(argv[5]);
        op          = atoi(argv[6]);

        // Ler dados da imagem .raw
        ReadImage rImage(filename, coluna, linha);
        st_image = rImage.vectorImage();
        st_image.vi_bits = atoi(argv[4]);

    } else {
        cout << "<file> <linha> <coluna> <Ng> <distancia>" << endl;
        //a.exit(EXIT_FAILURE);
        std::exit(EXIT_FAILURE);
    }

//    std::cout << "Dim: " << st_image.vi_linha << "x" << st_image.vi_coluna << endl ;
//    std::cout << "GL:   " << Ng << endl;
//    std::cout << endl;

    // **************************************** CPU ******************************************* //

    ///*
    double* matrizCoN_CPU = new double[Ng*Ng];
    Haralick haralick(&st_image);
    //MatrizCoocorrencia haralick(&st_image);

    CpuTimer Tempo;
//    startTime = omp_get_wtime();

    Tempo.Start();

    haralick.calcularMatrizCoN(matrizCoN_CPU, distancia); // double
    haralick.atCpu(matrizCoN_CPU, Ng);
    //haralick.CalcularAtributos();

    Tempo.Stop();
    Tempo.ElapsedSec();
//    stopTime = omp_get_wtime();

//    cout << "Cpu Time: " << stopTime - startTime << endl;

    //Plot Matriz
    //Img::ConstruirImgMC(matrizCoN_CPU, Ng);
    //cout << "**Gerada IMG da Matriz: " << endl;

    //*/

    // ************************************************************************************** //

//    atCpu atCpu(matrizCoN_CPU, Ng);

    double valor = 0.0;

    cout.unsetf ( std::ios::floatfield );
    cout.precision(15);
    cout << std::fixed << endl;

/*
    tempo::start();
    valor = haralick.max();
    tempo::stop();
    cout << "          MAX: " << valor << endl << endl;

    tempo::start();
    valor = haralick.min();
    tempo::stop();
    cout << "          MIN: " << valor << endl << endl;

    tempo::start();
    valor = haralick.media();
    tempo::stop();
    cout << "        Média: " << valor << endl << endl;

    tempo::start();
    valor = haralick.mediana();
    tempo::stop();
    cout << "      Mediana: " << valor << endl << endl;
    */

///*
    if (op == 1 || op == 0) {
//        cout << "Engergia" << endl;
        valor = haralick.energia();
        cout << "      Energia: " << valor << endl << endl;
    } if (op == 2 || op == 0) {
//        cout << "Contraste" << endl;
        valor = haralick.contraste();
        cout << "    Contraste: " << valor << endl << endl;
    } if (op == 3 || op == 0) {
//        cout << "Correlacao" << endl;
        valor = haralick.correlacao();
        cout << "   Correlacao: " << valor << endl << endl;
    } if (op == 4 || op == 0) {
//        cout << "Variancia" << endl;
        valor = haralick.variancia();
        cout << "    Variância: " << valor << endl << endl;
    } if (op == 5 || op == 0) {
//        cout << "MDI" << endl;
        valor = haralick.mdi();
        cout << "          MDI: " << valor << endl << endl;
    } if (op == 6 || op == 0) {
//        cout << "Entropia" << endl;
        valor = haralick.entropia();
        cout << "     Entropia: " << valor << endl << endl;
    } if (op == 7 || op == 0) {
//        cout << "Media Soma" << endl;
        valor = haralick.mediaSoma();
        cout << "   Media Soma: " << valor << endl << endl;
    } if (op == 8 || op == 0) {
//        cout << "Var Soma" << endl;
        valor = haralick.varianciaSoma();
        cout << "VarianciaSoma: " << valor << endl << endl;
    } if (op == 9 || op == 0) {
//        cout << "Soma Entr." << endl;
        valor = haralick.somaEntropia();
        cout << "Soma Entropia: " << valor << endl << endl;
    } if (op == 10 || op == 0) {
//        cout << "Var Dif." << endl;
        valor = haralick.varianciaDiferenca();
        cout << "Variancia Dif: " << valor << endl << endl;
    } if (op == 11 || op == 0) {
//        cout << "Dif Entr." << endl;
        valor = haralick.diferencaEntropia();
        cout << " Dif Entropia: " << valor << endl << endl;
    } if (op == 12 || op == 0) {
//        cout << "Med Correlacao 1" << endl;
        valor = haralick.medidasCorrelacao1();;
        cout << "Med Correlacao 1: " << valor << endl << endl;
    } if (op == 13 || op == 0){
//        cout << "Med Correlacao 2" << endl;
        valor = haralick.medidasCorrelacao2();
        cout << "Med Correlacao 2: " << valor << endl << endl;
    }
// */
    /*
    tempo::start();
    valor = haralick.mediaDiferenca();
    tempo::stop();
    cout << "    Media Dif: " << valor << endl << endl;
    */

    // */

    // ************************************************************************************** //

    //arq::arquivo("_cpu.txt", matrizCoN_cpu, nivelCinza);
    //arq::arquivo("_gpu.txt", matrizCoN_cpu, nivelCinza);

    // ************************************************************************************** //

    delete [] matrizCoN_CPU;

    // ************************************************************************************** //

    //return a.exec();
    //a.exit(EXIT_SUCCESS);
    return 0;
}
