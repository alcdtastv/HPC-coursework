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
        for (int row = 0; row < Ny; ++row)
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

            for (int col = 3; col < Nx - 3; ++col)
            {
                out[row + Ny * col] = -in[row + Ny * (col - 3)] / 60 + 3 * in[row + Ny * (col - 2)] / 20 - 3 * in[row + Ny * (col - 1)] / 4 +
                                      3 * in[row + Ny * (col + 1)] / 4 - 3 * in[row + Ny * (col + 2)] / 20 + in[row + Ny * (col + 3)] / 60;
            }
        }
    }
    void LoopDerivativeY(double *in, double *out)
    {
        for (int col = 0; col < Nx; ++col)
        {
            out[Ny * col] = -in[Ny * (col + 1) - 3] / 60 + 3 * in[Ny * (col + 1) - 2] / 20 - 3 * in[Ny * (col + 1) - 1] / 4 +
                            3 * in[Ny * col + 1] / 4 - 3 * in[Ny * col + 2] / 20 + in[Ny * col + 3] / 60;
            out[Ny * col + 1] = -in[Ny * (col + 1) - 2] / 60 + 3 * in[Ny * (col + 1) - 1] / 20 - 3 * in[Ny * col] / 4 +
                                3 * in[Ny * col + 2] / 4 - 3 * in[Ny * col + 3] / 20 + in[Ny * col + 4] / 60;
            out[Ny * col + 2] = -in[Ny * (col + 1) - 1] / 60 + 3 * in[Ny * col] / 20 - 3 * in[Ny * col + 1] / 4 +
                                3 * in[Ny * col + 3] / 4 - 3 * in[Ny * col + 4] / 20 + in[Ny * col + 5] / 60;
            out[Ny * (col + 1) - 3] = -in[Ny * (col + 1) - 6] / 60 + 3 * in[Ny * (col + 1) - 5] / 20 - 3 * in[Ny * (col + 1) - 4] / 4 +
                                      3 * in[Ny * (col + 1) - 2] / 4 - 3 * in[Ny * (col + 1) - 1] / 20 + in[Ny * col] / 60;
            out[Ny * (col + 1) - 2] = -in[Ny * (col + 1) - 5] / 60 + 3 * in[Ny * (col + 1) - 4] / 20 - 3 * in[Ny * (col + 1) - 3] / 4 +
                                      3 * in[Ny * (col + 1) - 1] / 4 - 3 * in[Ny * col] / 20 + in[Ny * col + 1] / 60;
            out[Ny * (col + 1) - 1] = -in[Ny * (col + 1) - 4] / 60 + 3 * in[Ny * (col + 1) - 3] / 20 - 3 * in[Ny * (col + 1) - 2] / 4 +
                                      3 * in[Ny * col] / 4 - 3 * in[Ny * col + 1] / 20 + in[Ny * col + 2] / 60;

            for (int row = 3; row < Ny - 3; ++row)
            {
                out[row + Ny * col] = -in[row - 3 + Ny * col] / 60 + 3 * in[row - 2 + Ny * col] / 20 - 3 * in[row - 1 + Ny * col] / 4 +
                                      3 * in[row + 1 + Ny * col] / 4 - 3 * in[row + 2 + Ny * col] / 20 + in[row + 3 + Ny * col] / 60;
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
    for (int col = 0; col < Nx; ++col)
        {
            out[Ny * col] = -in[Ny * (col + 1) - 3] / 60 + 3 * in[Ny * (col + 1) - 2] / 20 - 3 * in[Ny * (col + 1) - 1] / 4 +
                            3 * in[Ny * col + 1] / 4 - 3 * in[Ny * col + 2] / 20 + in[Ny * col + 3] / 60;
            out[Ny * col + 1] = -in[Ny * (col + 1) - 2] / 60 + 3 * in[Ny * (col + 1) - 1] / 20 - 3 * in[Ny * col] / 4 +
                                3 * in[Ny * col + 2] / 4 - 3 * in[Ny * col + 3] / 20 + in[Ny * col + 4] / 60;
            out[Ny * col + 2] = -in[Ny * (col + 1) - 1] / 60 + 3 * in[Ny * col] / 20 - 3 * in[Ny * col + 1] / 4 +
                                3 * in[Ny * col + 3] / 4 - 3 * in[Ny * col + 4] / 20 + in[Ny * col + 5] / 60;
            out[Ny * (col + 1) - 3] = -in[Ny * (col + 1) - 6] / 60 + 3 * in[Ny * (col + 1) - 5] / 20 - 3 * in[Ny * (col + 1) - 4] / 4 +
                                      3 * in[Ny * (col + 1) - 2] / 4 - 3 * in[Ny * (col + 1) - 1] / 20 + in[Ny * col] / 60;
            out[Ny * (col + 1) - 2] = -in[Ny * (col + 1) - 5] / 60 + 3 * in[Ny * (col + 1) - 4] / 20 - 3 * in[Ny * (col + 1) - 3] / 4 +
                                      3 * in[Ny * (col + 1) - 1] / 4 - 3 * in[Ny * col] / 20 + in[Ny * col + 1] / 60;
            out[Ny * (col + 1) - 1] = -in[Ny * (col + 1) - 4] / 60 + 3 * in[Ny * (col + 1) - 3] / 20 - 3 * in[Ny * (col + 1) - 2] / 4 +
                                      3 * in[Ny * col] / 4 - 3 * in[Ny * col + 1] / 20 + in[Ny * col + 2] / 60;

            for (int row = 3; row < Ny - 3; ++row)
            {
                out[row + Ny * col] = -in[row - 3 + Ny * col] / 60 + 3 * in[row - 2 + Ny * col] / 20 - 3 * in[row - 1 + Ny * col] / 4 +
                                      3 * in[row + 1 + Ny * col] / 4 - 3 * in[row + 2 + Ny * col] / 20 + in[row + 3 + Ny * col] / 60;
            }
        }
}

void printmatrix(double *matrix, int n1, int n2)
{
    for (int i = 0; i < n2; i++)
    {
        for (int j = 0; j < n1; j++)
        {
            cout << matrix[i + n2 * j] << ' ';
        }
        cout << endl;
    }
}

int main()
{
    double *in = new double[150];
    double *out = new double[150];

    for (int i = 0; i < 150; ++i)
    {
        in[i] = 1;
    }

    LoopDerivativeX(in, out, 15, 10);

    printmatrix(in, 15,10);
    printmatrix(out, 15,10);
}