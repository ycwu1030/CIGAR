#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;

double func_weight(double x, void* param)
{
    double dx = *((double *)param);
    double xmin = -1.0 + dx;
    double xmax = 1.0 - dx;
    double x_true = x*2.0-1.0;
    if (x_true < xmin || x_true > xmax)
    {
        return 0;
    }
    return (1.0+x_true*x_true)/(1.0-x_true*x_true)*2;
}
int main(int argc, char const *argv[])
{
    VEGAS_Integrator inter;
    
    for (double dx = 0.02; dx < 0.31; dx += 0.02)
    {
        inter.Set_Integrand(func_weight, &dx);
        inter.Improve_Grid();
        inter.Integration();
        cout<<"dx: "<<dx<<" res: "<<inter.Get_Result()<<" err: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<endl;
    }
    return 0;
}
