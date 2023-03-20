#include <iostream>
#include <fstream>
#include <iomanip>
#include <tgmath.h>
#include <boost/program_options.hpp>
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

class ShallowWater
{
public:
    void SetParameters(int argc, char **argv)
    {
        namespace po = boost::program_options;

        po::options_description desc("Allowed options");
        desc.add_options()("dt", po::value<double>(&dt), "Time-step to use.")("T", po::value<double>(&T), "Total integration time.")("Nx", po::value<int>(&Nx), "Number of grid points in x")("Ny", po::value<int>(&Ny), "Number of grid points in y")("ic", po::value<int>(&ic), "Index of the initial condition to use (1-4)")("type", po::value<char>(&type), "Blas or Loop");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        U = new double[Nx * Ny];
        V = new double[Nx * Ny];
        H = new double[Nx * Ny];

        if (type == 'B')
        {
            stencil = new double[7 * Nx]; // defined as x stencil, transpose if calculating y derivative

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx; ++i)
            {
                stencil[7 * i] = 1.0 / 60;
                stencil[7 * i + 1] = -3.0 / 20;
                stencil[7 * i + 2] = 3.0 / 4;
                stencil[7 * i + 3] = 0;
                stencil[7 * i + 4] = -3.0 / 4;
                stencil[7 * i + 5] = 3.0 / 20;
                stencil[7 * i + 6] = -1.0 / 60;
            }

            extrarowstencil = new double[Nx + 5];
            extrarowstencil[0] = 3.0 / 4;
            extrarowstencil[1] = -3.0 / 20;
            extrarowstencil[2] = 1.0 / 60;
            extrarowstencil[Nx + 2] = -1.0 / 60;
            extrarowstencil[Nx + 3] = 3.0 / 20;
            extrarowstencil[Nx + 4] = -3.0 / 4;

            extracolumnstencil = new double[Ny + 5];
            extracolumnstencil[0] = 3.0 / 4;
            extracolumnstencil[1] = -3.0 / 20;
            extracolumnstencil[2] = 1.0 / 60;
            extracolumnstencil[Ny + 2] = -1.0 / 60;
            extracolumnstencil[Ny + 3] = 3.0 / 20;
            extracolumnstencil[Ny + 4] = -3.0 / 4;
        }
    }
    void SetInitialConditions()
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
                    H[i + Ny * j] = 10 + exp(-(pow(i - 25, 2)+pow(j - 25, 2) )/ 25) + exp(-(pow(i - 75, 2)+pow(j - 75, 2) )/ 25);
                }
            }
        }
        else
        {
            cout << "Invalid initial condition" << endl;
        }
    }
    void xDerivative(double *in, double *out)
    {
        if (type == 'B')
        {
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Ny; ++i)
            {
                F77NAME(dgbmv)
                ('N', Nx, Nx, 3, 3, 1.0, stencil, 7, in + i, Nx, 0, out + i, Nx);

                out[i] += F77NAME(ddot)(Nx, extrarowstencil + 5, 1, in + i, Nx);
                out[i + Ny] += F77NAME(ddot)(Nx, extrarowstencil + 4, 1, in + i, Nx);
                out[i + 2 * Ny] += F77NAME(ddot)(Nx, extrarowstencil + 3, 1, in + i, Nx);
                out[i + Nx * (Ny - 3)] += F77NAME(ddot)(Nx, extrarowstencil + 2, 1, in + i, Nx);
                out[i + Nx * (Ny - 2)] += F77NAME(ddot)(Nx, extrarowstencil + 1, 1, in + i, Nx);
                out[i + Nx * (Ny - 1)] += F77NAME(ddot)(Nx, extrarowstencil, 1, in + i, Nx);
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

                for (int col = 3; col < Nx - 3; ++col)
                {
                    out[row + Ny * col] = -in[row + Ny * (col - 3)] * 0.0166667 + in[row + Ny * (col - 2)] * 0.15 - in[row + Ny * (col - 1)] * 0.75 +
                                          in[row + Ny * (col + 1)] * 0.75 - in[row + Ny * (col + 2)] *0.15 + in[row + Ny * (col + 3)] * 0.0166667;
                }
            }
        }
    }
    void yDerivative(double *in, double *out)
    {
        if (type == 'B')
        {
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Ny; ++i)
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
                                          
                for (int row = 3; row < Ny - 3; ++row)
                {
                    out[row + Ny * col] = -in[row - 3 + Ny * col] * 0.0166667 + in[row - 2 + Ny * col] * 0.15 - in[row - 1 + Ny * col] * 0.75 +
                                          in[row + 1 + Ny * col] * 0.75 - in[row + 2 + Ny * col] * 0.15 + in[row + 3 + Ny * col] * 0.0166667;
                }
            }
        }
    }
    void TimeIntegrate()
    {
        // Initialising variables
        double *u = new double[Nx * Ny];
        double *v = new double[Nx * Ny];
        double *h = new double[Nx * Ny];
        double *hu = new double[Nx * Ny];
        double *hv = new double[Nx * Ny];

        double *dudx = new double[Nx * Ny];
        double *dudy = new double[Nx * Ny];

        double *dvdx = new double[Nx * Ny];
        double *dvdy = new double[Nx * Ny];

        double *dhdx = new double[Nx * Ny];
        double *dhdy = new double[Nx * Ny];

        double *dhudx = new double[Nx * Ny];
        double *dhvdy = new double[Nx * Ny];

        double *k1u = new double[Nx * Ny];
        double *k1v = new double[Nx * Ny];
        double *k1h = new double[Nx * Ny];

        double *k2u = new double[Nx * Ny];
        double *k2v = new double[Nx * Ny];
        double *k2h = new double[Nx * Ny];

        double *k3u = new double[Nx * Ny];
        double *k3v = new double[Nx * Ny];
        double *k3h = new double[Nx * Ny];

        double *k4u = new double[Nx * Ny];
        double *k4v = new double[Nx * Ny];
        double *k4h = new double[Nx * Ny];

        double g = 9.81;

        for (int j = 0; j < T / dt; ++j)
        {
            // Copying U, V, H to u, v, h (temp arrays) and calculating hu, hv
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i];
                v[i] = V[i];
                h[i] = H[i];
                hu[i] = H[i] * U[i];
                hv[i] = H[i] * V[i];
            }

            // Calculating K1
            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);
            xDerivative(hu, dhudx);
            yDerivative(hv, dhvdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k1u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k1v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k1h[i] = -dhudx[i] - dhvdy[i];
            }

            // Calculating K2
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k1u[i] * dt / 2;
                v[i] = V[i] + k1v[i] * dt / 2;
                h[i] = H[i] + k1h[i] * dt / 2;
                hu[i] = h[i] * u[i];
                hv[i] = h[i] * v[i];
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);
            xDerivative(hu, dhudx);
            yDerivative(hv, dhvdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k2u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k2v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k2h[i] = -dhudx[i] - dhvdy[i];
            }

            // Calculating K3
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k2u[i] * dt / 2;
                v[i] = V[i] + k2v[i] * dt / 2;
                h[i] = H[i] + k2h[i] * dt / 2;
                hu[i] = h[i] * u[i];
                hv[i] = h[i] * v[i];
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);
            xDerivative(hu, dhudx);
            yDerivative(hv, dhvdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k3u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k3v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k3h[i] = -dhudx[i] - dhvdy[i];
            }

            // Calculating K4
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                u[i] = U[i] + k3u[i] * dt;
                v[i] = V[i] + k3v[i] * dt;
                h[i] = H[i] + k3h[i] * dt;
                hu[i] = h[i] * u[i];
                hv[i] = h[i] * v[i];
            }

            xDerivative(u, dudx);
            yDerivative(u, dudy);
            xDerivative(v, dvdx);
            yDerivative(v, dvdy);
            xDerivative(h, dhdx);
            yDerivative(h, dhdy);
            xDerivative(hu, dhudx);
            yDerivative(hv, dhvdy);

            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                k4u[i] = -u[i] * dudx[i] - v[i] * dudy[i] - g * dhdx[i];
                k4v[i] = -u[i] * dvdx[i] - v[i] * dvdy[i] - g * dhdy[i];
                k4h[i] = -dhudx[i] - dhvdy[i];
            }

            // Calculating U, V, and H
            #pragma omp parallel for default(shared) schedule(static)
            for (int i = 0; i < Nx * Ny; ++i)
            {
                U[i] = U[i] + (k1u[i] + 2 * k2u[i] + 2 * k3u[i] + k4u[i]) * dt / 6;
                V[i] = V[i] + (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]) * dt / 6;
                H[i] = H[i] + (k1h[i] + 2 * k2h[i] + 2 * k3h[i] + k4h[i]) * dt / 6;
            }
        }
    }
    void output()
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

private:
    double dt;
    double T;
    int Nx;
    int Ny;
    int ic;
    char type;
    double *U;
    double *V;
    double *H;
    double *stencil;
    double *extrarowstencil;
    double *extracolumnstencil;
};

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