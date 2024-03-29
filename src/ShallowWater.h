/**
 * @file ShallowWater.h
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the class definition of the ShallowWater class.
 * @version 1.0
 * @date 2023-03-20
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef ShallowWater_H
#define ShallowWater_H

using namespace std;

/**
 * @class ShallowWater
 * @brief Defines variables and methods to implement a shallow water equation solver.
 *
 */

class ShallowWater
{
public:
    ShallowWater(int argc, char **argv);
    
    ~ShallowWater();

    void SetInitialConditions();

    void xDerivative(double *in, double *out);

    void yDerivative(double *in, double *out);

    void TimeIntegrate();

    void output();

private:
    double dt;                  ///< Time step
    double T;                   ///< Total integration time
    int Nx;                     ///< Number of x grid points
    int Ny;                     ///< Number of y grid points
    int ic;                     ///< Initial condition selector
    char type;                  ///< Loop or Blas selector
    double *U;                  ///< X component of velocity
    double *V;                  ///< Y component of velocity
    double *H;                  ///< Height
    double *stencil;            ///< Stencil for blas derivative calculation
    double *extrarowstencil;    ///< Stencil for extra rows not accounted by the dgbmv operation
    double *extracolumnstencil; ///< Stencil for extra columns not accounted by the dgbmv operation
};

#endif