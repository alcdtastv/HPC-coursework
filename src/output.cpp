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

    for (int row = 0; row < Ny; ++row)
    {
        for (int col = 0; col < Nx; ++col)
        {
            out << col << ' ' << row << ' ' << U[row + Ny * col] << ' ' << V[row + Ny * col] << ' ' << H[row + Ny * col] << endl;
        }
        out << endl;
    }

    out.close();
}