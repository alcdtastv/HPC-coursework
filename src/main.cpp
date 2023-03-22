/**
 * @file main.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This program solves the shallow water equation using runge kutta and a central differencing scheme,
 *        implemented through both a loop * based method and blas.
 * @version 0.1
 * @date 2023-03-20
 *
 * @copyright Copyright (c) 2023
 *
 */

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