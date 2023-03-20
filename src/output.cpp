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