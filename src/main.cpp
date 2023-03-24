/**
 * @file main.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the main function of the shallow water equations solver.
 * @version 1.0
 * @date 2023-03-20
 *
 * @copyright Copyright (c) 2023
 *
 */
//l
#include "ShallowWater.h"
#include <omp.h>
#include <iostream>

using namespace std;

/**
 * @brief 
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */

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

    ShallowWater wave(argc, argv);

    wave.SetInitialConditions();

    wave.TimeIntegrate();

    wave.output();
}