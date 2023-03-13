#include <iostream>

using namespace std;

class ShallowWater
{
public:
    void SetInitialConditions()
    {
        return;
    }
    void SetParameters()
    {
        return;
    }
    void TimeIntegrate()
    {
        return;
    }
    void BlasDerivativeX()
    {
        return;
    }
    void LoopDerivativeX(double *in, double *out)
    {
        for (int row = 0; row < Nx; ++row)
        {
            out[row] = -in[(Nx - 3) * Ny + row] / 60 + 3 * in[(Nx - 2) * Ny + row] / 20 - 3 * in[(Nx - 1) * Ny + row] / 4 +
                       3 * in[Ny + row] / 4 - 3 * in[2 * Ny + row] / 20 + in[3 * Ny + row] / 60;

            out[Ny + row] = -in[(Nx - 2) * Ny + row] / 60 + 3 * in[(Nx - 1) * Ny + row] / 20 - 3 * in[row] / 4 +
                            3 * in[2 * Ny + row] / 4 - 3 * in[3 * Ny + row] / 20 + in[4 * Ny + row] / 60;

            out[2 * Ny + row] = -in[(Nx - 1) * Ny + row] / 60 + 3 * in[row] / 20 - 3 * in[Ny + row] / 4 +
                                3 * in[3 * Ny + row] / 4 - 3 * in[4 * Ny + row] / 20 + in[5 * Ny + row] / 60;

            out[(Nx - 3) * Ny + row] = -in[(Nx - 6) * Ny + row] / 60 + 3 * in[(Nx - 5) * Ny + row] / 20 - 3 * in[(Nx - 4) * Ny + row] / 4 +
                                       3 * in[(Nx - 2) * Ny + row] / 4 - 3 * in[(Nx - 1) * Ny + row] / 20 + in[row] / 60;

            out[(Nx - 2) * Ny + row] = -in[(Nx - 5) * Ny + row] / 60 + 3 * in[(Nx - 4) * Ny + row] / 20 - 3 * in[(Nx - 3) * Ny + row] / 4 +
                                       3 * in[(Nx - 1) * Ny + row] / 4 - 3 * in[row] / 20 + in[Ny + row] / 60;

            out[(Nx - 1) * Ny + row] = -in[(Nx - 4) * Ny + row] / 60 + 3 * in[(Nx - 3) * Ny + row] / 20 - 3 * in[(Nx - 2) * Ny + row] / 4 +
                                       3 * in[row] / 4 - 3 * in[Ny + row] / 20 + in[2 * Ny + row] / 60;

            for (int col = 3; col < Ny - 3; ++col)
            {
                out[row + Ny * col] = -in[row + Ny * (col - 3)] / 60 + 3 * in[row + Ny * (col - 2)] / 20 - 3 * in[row + Ny * (col - 1)] / 4 +
                                      3 * in[row + Ny * (col + 1)] / 4 - 3 * in[row + Ny * (col + 2)] / 20 + in[row + Ny * (col + 3)] / 60;
            }
        }
    }

private:
    double dt = 0.001;
    double T = 10;
    int Nx = 100;
    int Ny = 100;
    int ic;
    double *U;
    double *V;
    double *H;
};

void LoopDerivativeX(double *in, double *out, int Nx, int Ny)
{
    for (int row = 0; row < Nx; ++row)
    {
        out[row] = -in[(Nx - 3) * Ny + row] / 60 + 3 * in[(Nx - 2) * Ny + row] / 20 - 3 * in[(Nx - 1) * Ny + row] / 4 +
                   3 * in[Ny + row] / 4 - 3 * in[2 * Ny + row] / 20 + in[3 * Ny + row] / 60;

        out[Ny + row] = -in[(Nx - 2) * Ny + row] / 60 + 3 * in[(Nx - 1) * Ny + row] / 20 - 3 * in[row] / 4 +
                        3 * in[2 * Ny + row] / 4 - 3 * in[3 * Ny + row] / 20 + in[4 * Ny + row] / 60;

        out[2 * Ny + row] = -in[(Nx - 1) * Ny + row] / 60 + 3 * in[row] / 20 - 3 * in[Ny + row] / 4 +
                            3 * in[3 * Ny + row] / 4 - 3 * in[4 * Ny + row] / 20 + in[5 * Ny + row] / 60;

        out[(Nx - 3) * Ny + row] = -in[(Nx - 6) * Ny + row] / 60 + 3 * in[(Nx - 5) * Ny + row] / 20 - 3 * in[(Nx - 4) * Ny + row] / 4 +
                                   3 * in[(Nx - 2) * Ny + row] / 4 - 3 * in[(Nx - 1) * Ny + row] / 20 + in[row] / 60;

        out[(Nx - 2) * Ny + row] = -in[(Nx - 5) * Ny + row] / 60 + 3 * in[(Nx - 4) * Ny + row] / 20 - 3 * in[(Nx - 3) * Ny + row] / 4 +
                                   3 * in[(Nx - 1) * Ny + row] / 4 - 3 * in[row] / 20 + in[Ny + row] / 60;

        out[(Nx - 1) * Ny + row] = -in[(Nx - 4) * Ny + row] / 60 + 3 * in[(Nx - 3) * Ny + row] / 20 - 3 * in[(Nx - 2) * Ny + row] / 4 +
                                   3 * in[row] / 4 - 3 * in[Ny + row] / 20 + in[2 * Ny + row] / 60;

        for (int col = 3; col < Ny - 3; ++col)
        {
            out[row + Ny * col] = -in[row + Ny * (col - 3)] / 60 + 3 * in[row + Ny * (col - 2)] / 20 - 3 * in[row + Ny * (col - 1)] / 4 +
                                  3 * in[row + Ny * (col + 1)] / 4 - 3 * in[row + Ny * (col + 2)] / 20 + in[row + Ny * (col + 3)] / 60;
        }
    }
}

void printmatrix(double *matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << matrix[i + n * j] << ' ';
        }
        cout << endl;
    }
}

int main()
{
    double *in = new double[100];
    double *out = new double[100];

    for (int i = 0; i < 100; ++i)
    {
        in[i] = 1;
    }

    LoopDerivativeX(in, out, 10, 10);

    printmatrix(in, 10);
    printmatrix(out, 10);
}