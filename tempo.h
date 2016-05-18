#ifndef TEMPO
#define TEMPO

#include <stdio.h>
//#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

struct CpuTimer
{
    double start;
    double stop;

    CpuTimer() {}

    ~CpuTimer() {}

    inline void Start()
    {
        start = omp_get_wtime();
    }

    inline void Stop()
    {
        stop = omp_get_wtime();
    }

    inline double ElapsedSec()
    {
        double tempo = stop - start;
        printf("%f\n", tempo);

        return tempo;
    }
};
/*
struct timeval timerStart;

void StartTimer()
{
    gettimeofday(&timerStart, NULL);
}

// time elapsed in ms
double GetTimer()
{
    struct timeval timerStop, timerElapsed;
    gettimeofday(&timerStop, NULL);
    timersub(&timerStop, &timerStart, &timerElapsed);
    return timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
}
*/
#endif // TEMPO

