#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;

double func_weight(vector<double> x, void* param)
{
    return 0;
}
int main(int argc, char const *argv[])
{
    VEGAS_Integrator inter;
    
    double dx = 0.2;
    inter.Set_Integrand(func_weight, 1, &dx);
    inter.Improve_Grid();
    inter.Integration();
    cout<<"dx: "<<dx<<" res: "<<inter.Get_Result()<<" err: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<endl;
    return 0;
}
