#ifndef UTILIDADES
#define UTILIDADES

#include <stdio.h>
#include <stdlib.h>

using namespace std;

namespace cmp {

    bool compara(double *mCpu, double *mGpu, int tamM) {

        std::cout << "cpu != gpu" << endl;
        std::cout.unsetf ( std::ios::floatfield );
        std::cout.precision(15);

        for (int i = 0; i < tamM; i++) {
            for (int j = 0; j < tamM; j++) {
                if (mCpu[i * tamM + j] != mGpu[i * tamM + j]) {
                    //printf("%d ,%d : %f != %lf\n", i, j, mCpu[i * tamM + j], mGpu[i * tamM + j]);
                    /*
                    std::cout << i << ", " << j << ": " << std::fixed
                              << mCpu[i * tamM + j] << " != "
                              << mGpu[i * tamM + j] << endl;
                    */
                    double result = mCpu[i * tamM + j] - mGpu[i * tamM + j];
                    if (result < 0.0) {
                        result = result * (-1.0);
                    }
                    if ( result > 0.0000000001)
                        std::cout << std::fixed << result << std::endl;
                    //return false;
                    //if ( (i*tamM + j) == 1024) return false;
                }
            }
        }

        return true;
    }
}

namespace Img {
    void ConstruirImgMC(double *matriz, int TAM) {

        char filename[] = "matriz.pgm";

        remove(filename);

        FILE *file = fopen(filename, "w");
        if (file == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }

        fprintf(file, "P2\n");
        fprintf(file, "%d %d\n", TAM, TAM);
        fprintf(file, "255\n");
        for(int i = 0; i < TAM*TAM; ++i) {
            double nv  = matriz[i];
            if (nv > 0.0)
                fprintf(file, "0\n");
            else
                fprintf(file, "255\n");
        }

        fclose(file);
    }
}

namespace tempo {

    double startTime = 0.0;
    double stopTime = 0.0;

    inline void start() {
        startTime = omp_get_wtime();
    }

    inline void stop() {
        stopTime = omp_get_wtime();
        std::cout << "         Time: " << stopTime - startTime << std::endl;
    }
}
#endif // UTILIDADES

