#include <iostream>
#include <tgmath.h>
#include <boost/program_options.hpp>

using namespace std;

class ShallowWater
{
public:
    void SetParameters(int argc, char **argv)
    {
        namespace po = boost::program_options;

        po::options_description desc("Allowed options");
        desc.add_options()("dt", po::value<double>(&dt), "Time-step to use.")("T", po::value<double>(&T), "Total integration time.")("Nx", po::value<int>(&Nx), "Number of grid points in x")("Ny", po::value<int>(&Ny), "Number of grid points in y")("ic", po::value<int>(&ic), "Index of the initial condition to use (1-4)");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        U = new double[Nx * Ny];
        V = new double[Nx * Ny];
        H = new double[Nx * Ny];
    }
    void SetInitialConditions()
    {
        if (ic == 1)
        {
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
            for (int i = 0; i < Nx; ++i)
            {
                for (int j = 0; j < Ny; ++j)
                {
                    H[i + Ny * j] = 10 + exp(-pow(i - 25, 2) / 25) + exp(-pow(j - 25, 2) / 25) + exp(-pow(i - 75, 2) / 25) + exp(-pow(j - 75, 2) / 25);
                }
            }
        }
        else
        {
            cout << "Invalid initial condition" << endl;
        }
    }
    void BlasDerivativeX()
    {
        return;
    }
    void BlasDerivativeY()
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
    void TimeIntegrate()
    {
        
    }
    void print()
    {
        cout << "dt = " << dt << endl;
        cout << "T = " << T << endl;
        cout << "Nx = " << Nx << endl;
        cout << "Ny = " << Ny << endl;
        cout << "ic = " << ic << endl;
    }

    void printmatrix(double *matrix)
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                cout << matrix[i + Ny * j] << ' ';
            }
            cout << endl;
        }
    }

    double *getH()
    {
        return H;
    }

private:
    double dt;
    double T;
    int Nx;
    int Ny;
    int ic;
    double *U;
    double *V;
    double *H;
};



int main(int argc, char **argv)
{
    ShallowWater obj;

    obj.SetParameters(argc, argv);

    obj.print();

    obj.SetInitialConditions();

    obj.printmatrix(obj.getH());
}