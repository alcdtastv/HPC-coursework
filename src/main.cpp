#include "ShallowWater.h"
#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    int nthreads, threadid;
    #pragma omp parallel private(threadid)
    {
        threadid = omp_get_thread_num();
        if (threadid == 0)
        {
            nthreads = omp_get_num_threads();
            cout << "Number of threads: " << nthreads << endl;
        }
    }
    
    ShallowWater wave;

    wave.SetParameters(argc, argv);

    wave.SetInitialConditions();

    wave.TimeIntegrate();

    wave.output();
}