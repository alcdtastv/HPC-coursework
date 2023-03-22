/**
 * @file SetInitialConditions.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-22
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "ShallowWater.h"
#include <omp.h>
#include <tgmath.h>
#include <iostream>

using namespace std;
/**
 * @brief Sets initial conditions for the H array based on the ic class variable.
 *
 */
void ShallowWater::SetInitialConditions()
{
    if (ic == 1)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                H[i + Ny * j] = 10 + exp(-pow(i - 50, 2) / 25);
            }
        }
    }
    else if (ic == 2)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                H[i + Ny * j] = 10 + exp(-pow(j - 50, 2) / 25);
            }
        }
    }
    else if (ic == 3)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                H[i + Ny * j] = 10 + exp(-(pow(i - 50, 2) + pow(j - 50, 2)) / 25);
            }
        }
    }
    else if (ic == 4)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Nx; ++i)
        {
            for (int j = 0; j < Ny; ++j)
            {
                H[i + Ny * j] = 10 + exp(-(pow(i - 25, 2) + pow(j - 25, 2)) / 25) + exp(-(pow(i - 75, 2) + pow(j - 75, 2)) / 25);
            }
        }
    }
    else
    {
        cout << "Invalid initial condition" << endl;
    }
}