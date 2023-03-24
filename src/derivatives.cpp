/**
 * @file derivatives.cpp
 * @author Luca Mazzotta (luca.mazzotta19@imperial.ac.uk)
 * @brief This file contains the functions to calculate the x and y derivatives.
 * @version 1.0
 * @date 2023-03-20
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "ShallowWater.h"
#include <omp.h>

#define F77NAME(x) x##_

using namespace std;

extern "C"
{
    void F77NAME(dgbmv)(const char &trans, const int &m, const int &n, const int &KL, const int &KU, const double &alpha,
                        const double *A, const int &lda, const double *x, const int &incx, const double &beta, double *y,
                        const int &incy);

    double F77NAME(ddot)(const int &n, const double *x, const int &incx, const double *y, const int &incy);
}

/**
 * @brief Calculates the x derivative of the input array using a central differencing scheme, implemented both through blas
 *        and loops (chosen based on the type class variable).
 *
 * @param in    Input array
 * @param out   Output array
 */

void ShallowWater::xDerivative(double *in, double *out)
{
    if (type == 'B')
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Ny; ++i)
        {
            F77NAME(dgbmv)
            ('N', Nx, Nx, 3, 3, 1.0, stencil, 7, in + i, Ny, 0, out + i, Ny);

            out[i] += F77NAME(ddot)(Nx, extrarowstencil + 5, 1, in + i, Ny);
            out[i + Ny] += F77NAME(ddot)(Nx, extrarowstencil + 4, 1, in + i, Ny);
            out[i + 2 * Ny] += F77NAME(ddot)(Nx, extrarowstencil + 3, 1, in + i, Ny);
            out[i + Ny * (Nx - 3)] += F77NAME(ddot)(Nx, extrarowstencil + 2, 1, in + i, Ny);
            out[i + Ny * (Nx - 2)] += F77NAME(ddot)(Nx, extrarowstencil + 1, 1, in + i, Ny);
            out[i + Ny * (Nx - 1)] += F77NAME(ddot)(Nx, extrarowstencil, 1, in + i, Ny);
        }
    }
    else if (type == 'L')
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int row = 0; row < Ny; ++row)
        {
            out[row] = -in[(Nx - 3) * Ny + row] * 0.0166667 + in[(Nx - 2) * Ny + row] * 0.15 - in[(Nx - 1) * Ny + row] * 0.75 +
                       in[Ny + row] * 0.75 - in[2 * Ny + row] * 0.15 + in[3 * Ny + row] * 0.0166667;

            out[Ny + row] = -in[(Nx - 2) * Ny + row] * 0.0166667 + in[(Nx - 1) * Ny + row] * 0.15 - in[row] * 0.75 +
                            in[2 * Ny + row] * 0.75 - in[3 * Ny + row] * 0.15 + in[4 * Ny + row] * 0.0166667;

            out[2 * Ny + row] = -in[(Nx - 1) * Ny + row] * 0.0166667 + in[row] * 0.15 - in[Ny + row] * 0.75 +
                                in[3 * Ny + row] * 0.75 - in[4 * Ny + row] * 0.15 + in[5 * Ny + row] * 0.0166667;

            out[(Nx - 3) * Ny + row] = -in[(Nx - 6) * Ny + row] * 0.0166667 + in[(Nx - 5) * Ny + row] * 0.15 - in[(Nx - 4) * Ny + row] * 0.75 +
                                       in[(Nx - 2) * Ny + row] * 0.75 - in[(Nx - 1) * Ny + row] * 0.15 + in[row] * 0.0166667;

            out[(Nx - 2) * Ny + row] = -in[(Nx - 5) * Ny + row] * 0.0166667 + in[(Nx - 4) * Ny + row] * 0.15 - in[(Nx - 3) * Ny + row] * 0.75 +
                                       in[(Nx - 1) * Ny + row] * 0.75 - in[row] * 0.15 + in[Ny + row] * 0.0166667;

            out[(Nx - 1) * Ny + row] = -in[(Nx - 4) * Ny + row] * 0.0166667 + in[(Nx - 3) * Ny + row] * 0.15 - in[(Nx - 2) * Ny + row] * 0.75 +
                                       in[row] * 0.75 - in[Ny + row] * 0.15 + in[2 * Ny + row] * 0.0166667;

            #pragma omp simd
            for (int col = 3; col < Nx - 3; ++col)
            {
                out[row + Ny * col] = -in[row + Ny * (col - 3)] * 0.0166667 + in[row + Ny * (col - 2)] * 0.15 - in[row + Ny * (col - 1)] * 0.75 +
                                      in[row + Ny * (col + 1)] * 0.75 - in[row + Ny * (col + 2)] * 0.15 + in[row + Ny * (col + 3)] * 0.0166667;
            }
        }
    }
}

/**
 * @brief Calculates the y derivative of the input array using a central differencing scheme, implemented both through blas
 *        and loops (chosen based on the type class variable).
 * 
 * @param in Input array
 * @param out Output array
 */

void ShallowWater::yDerivative(double *in, double *out)
{
    if (type == 'B')
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int i = 0; i < Nx; ++i)
        {
            F77NAME(dgbmv)
            ('T', Ny, Ny, 3, 3, -1.0, stencil, 7, in + i * Ny, 1, 0, out + i * Ny, 1);

            out[Ny * i] += F77NAME(ddot)(Ny, extracolumnstencil + 5, 1, in + i * Ny, 1);
            out[Ny * i + 1] += F77NAME(ddot)(Ny, extracolumnstencil + 4, 1, in + i * Ny, 1);
            out[Ny * i + 2] += F77NAME(ddot)(Ny, extracolumnstencil + 3, 1, in + i * Ny, 1);
            out[Ny * (i + 1) - 3] += F77NAME(ddot)(Ny, extracolumnstencil + 2, 1, in + i * Ny, 1);
            out[Ny * (i + 1) - 2] += F77NAME(ddot)(Ny, extracolumnstencil + 1, 1, in + i * Ny, 1);
            out[Ny * (i + 1) - 1] += F77NAME(ddot)(Ny, extracolumnstencil, 1, in + i * Ny, 1);
        }
    }
    else if (type == 'L')
    {
        #pragma omp parallel for default(shared) schedule(static)
        for (int col = 0; col < Nx; ++col)
        {
            out[Ny * col] = -in[Ny * (col + 1) - 3] * 0.0166667 + in[Ny * (col + 1) - 2] * 0.15 - in[Ny * (col + 1) - 1] * 0.75 +
                            in[Ny * col + 1] * 0.75 - in[Ny * col + 2] * 0.15 + in[Ny * col + 3] * 0.0166667;

            out[Ny * col + 1] = -in[Ny * (col + 1) - 2] * 0.0166667 + in[Ny * (col + 1) - 1] * 0.15 - in[Ny * col] * 0.75 +
                                in[Ny * col + 2] * 0.75 - in[Ny * col + 3] * 0.15 + in[Ny * col + 4] * 0.0166667;

            out[Ny * col + 2] = -in[Ny * (col + 1) - 1] * 0.0166667 + in[Ny * col] * 0.15 - in[Ny * col + 1] * 0.75 +
                                in[Ny * col + 3] * 0.75 - in[Ny * col + 4] * 0.15 + in[Ny * col + 5] * 0.0166667;

            out[Ny * (col + 1) - 3] = -in[Ny * (col + 1) - 6] * 0.0166667 + in[Ny * (col + 1) - 5] * 0.15 - in[Ny * (col + 1) - 4] * 0.75 +
                                      in[Ny * (col + 1) - 2] * 0.75 - in[Ny * (col + 1) - 1] * 0.15 + in[Ny * col] * 0.0166667;

            out[Ny * (col + 1) - 2] = -in[Ny * (col + 1) - 5] * 0.0166667 + in[Ny * (col + 1) - 4] * 0.15 - in[Ny * (col + 1) - 3] * 0.75 +
                                      in[Ny * (col + 1) - 1] * 0.75 - in[Ny * col] * 0.15 + in[Ny * col + 1] * 0.0166667;

            out[Ny * (col + 1) - 1] = -in[Ny * (col + 1) - 4] * 0.0166667 + in[Ny * (col + 1) - 3] * 0.15 - in[Ny * (col + 1) - 2] * 0.75 +
                                      in[Ny * col] * 0.75 - in[Ny * col + 1] * 0.15 + in[Ny * col + 2] * 0.0166667;

            #pragma omp simd
            for (int row = 3; row < Ny - 3; ++row)
            {
                out[row + Ny * col] = -in[row - 3 + Ny * col] * 0.0166667 + in[row - 2 + Ny * col] * 0.15 - in[row - 1 + Ny * col] * 0.75 +
                                      in[row + 1 + Ny * col] * 0.75 - in[row + 2 + Ny * col] * 0.15 + in[row + 3 + Ny * col] * 0.0166667;
            }
        }
    }
}