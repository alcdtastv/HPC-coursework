/**
 * @file output.cpp
 * @author Luca Mazzotta
 * @brief Outputs the U, V and H values to output.txt, just afer their x and y coordinates.
 * @version 0.1
 * @date 2023-03-22
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "ShallowWater.h"
#include <fstream>

using namespace std;

void ShallowWater::output()
{
    ofstream out("output.txt", ios::out | ios::trunc);

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            out << i << ' ' << j << ' ' << U[i + Ny * j] << ' ' << V[i + Ny * j] << ' ' << H[i + Ny * j] << endl;
        }
        out << endl;
    }

    out.close();
}