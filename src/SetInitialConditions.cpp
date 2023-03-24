/**
 * @file SetInitialConditions.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the function to set the initial conditions for the H array.
 * @version 1.0
 * @date 2023-03-20
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
        for (int row = 0; row < Ny; ++row)
        {
            for (int col = 0; col < Nx; ++col)
            {
                H[row + Ny * col] = 10 + exp(-pow(col - 50, 2) / 25);
            }
        }
    }
    else if (ic == 2)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int row = 0; row < Ny; ++row)
        {
            for (int col = 0; col < Nx; ++col)
            {
                H[row + Ny * col] = 10 + exp(-pow(row - 50, 2) / 25);
            }
        }
    }
    else if (ic == 3)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int row = 0; row < Ny; ++row)
        {
            for (int col = 0; col < Nx; ++col)
            {
                H[row + Ny * col] = 10 + exp(-(pow(row - 50, 2) + pow(col - 50, 2)) / 25);
            }
        }
    }
    else if (ic == 4)
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int row = 0; row < Ny; ++row)
        {
            for (int col = 0; col < Nx; ++col)
            {
                H[row + Ny * col] = 10 + exp(-(pow(col - 25, 2) + pow(row - 25, 2)) / 25) + exp(-(pow(col - 75, 2) + pow(row - 75, 2)) / 25);
            }
        }
    }
    else
    {
        cout << "Invalid initial condition" << endl;
    }
}