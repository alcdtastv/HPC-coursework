#ifndef ShallowWater_H
#define ShallowWater_H

using namespace std;

class ShallowWater
{
public:
    void SetParameters(int argc, char **argv);
    
    void SetInitialConditions();

    void xDerivative(double *in, double *out);
    
    void yDerivative(double *in, double *out);
    
    void TimeIntegrate();

    void output();

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

#endif