/**
 * @file ConstructorDestructor.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the constructor and destructor of the ShallowWater class.
 * @version 1.0
 * @date 2023-03-20
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ShallowWater.h"
#include <boost/program_options.hpp>

using namespace std;

/**
 * @brief Constructs the shallow water object and sets its parameters.
 *
 * @param argc
 * @param argv
 */

ShallowWater::ShallowWater(int argc, char **argv)
{
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("dt", po::value<double>(&dt), "Time-step to use.")
    ("T", po::value<double>(&T), "Total integration time.")
    ("Nx", po::value<int>(&Nx), "Number of grid points in x")
    ("Ny", po::value<int>(&Ny), "Number of grid points in y")
    ("ic", po::value<int>(&ic), "Index of the initial condition to use (1-4)")
    ("type", po::value<char>(&type), "Blas or Loop");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    U = new double[Nx * Ny];
    V = new double[Nx * Ny];
    H = new double[Nx * Ny];

    if (type == 'B')
    {
        int temp;
        if (Nx > Ny)
        {
            temp = Nx;
        }
        else
        {
            temp = Ny;
        }

        stencil = new double[7 * temp]; // sized to the largest domain dimension

        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < temp; ++i)
        {
            stencil[7 * i] = 1.0 / 60;
            stencil[7 * i + 1] = -3.0 / 20;
            stencil[7 * i + 2] = 3.0 / 4;
            stencil[7 * i + 3] = 0;
            stencil[7 * i + 4] = -3.0 / 4;
            stencil[7 * i + 5] = 3.0 / 20;
            stencil[7 * i + 6] = -1.0 / 60;
        }

        // stencil to account for the row which cannot be computed by dgbmv
        extrarowstencil = new double[Nx + 5];
        extrarowstencil[0] = 3.0 / 4;
        extrarowstencil[1] = -3.0 / 20;
        extrarowstencil[2] = 1.0 / 60;
        extrarowstencil[Nx + 2] = -1.0 / 60;
        extrarowstencil[Nx + 3] = 3.0 / 20;
        extrarowstencil[Nx + 4] = -3.0 / 4;

        // stencil to account for the column which cannot be computed by dgbmv
        extracolumnstencil = new double[Ny + 5];
        extracolumnstencil[0] = 3.0 / 4;
        extracolumnstencil[1] = -3.0 / 20;
        extracolumnstencil[2] = 1.0 / 60;
        extracolumnstencil[Ny + 2] = -1.0 / 60;
        extracolumnstencil[Ny + 3] = 3.0 / 20;
        extracolumnstencil[Ny + 4] = -3.0 / 4;
    }
}

ShallowWater::~ShallowWater()
{
    delete[] U;
    delete[] V;
    delete[] H;

    if (type == 'B')
    {
        delete[] stencil;
        delete[] extrarowstencil;
        delete[] extracolumnstencil;
    }
}