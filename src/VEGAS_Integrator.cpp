#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;


void VEGAS_Integrator::Set_Integrand(INTEGRAND integrand, void* param)
{
    func = integrand;
    userdata = param;
    Results.clear();
    Sigma2.clear();
    map.Reset_Map();
}

void VEGAS_Integrator::Improve_Grid(int Iter, int Neval)
{
    double y; // Random number between 0 to 1;
    double x; // The argument for integrand;
    double f_eval; // evaluated integrand value;
    for (int it = 0; it < Iter; it++)
    {
        for (int ne = 0; ne < Neval; ne++)
        {
            y = dist(rng);
            x = map.Get_X(y);
            f_eval = func(x,userdata);
            map.Accumulate_Weight(y,f_eval);  
        }
        map.Update_Map();
    }
}
void VEGAS_Integrator::Integration(int Iter, int Neval)
{
    double y; // Random number between 0 to 1;
    double x; // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; // The Jacobian from y to x;
    for (int it = 0; it < Iter; it++)
    {
        double res = 0;
        double sig2 = 0;
        for (int ne = 0; ne < Neval; ne++)
        {
            y = dist(rng);
            x = map.Get_X(y);
            f_eval = func(x,userdata);
            Jac = map.Get_Jac(y);
            res += f_eval*Jac;
            sig2 += pow(f_eval*Jac,2);
        }
        Results.push_back(res/Neval);
        Sigma2.push_back((sig2/Neval-pow(res/Neval,2))/(Neval-1.0));
    }
}
double VEGAS_Integrator::Get_Result()
{
    double res_num = 0;
    double res_den = 0;
    for (int i = 0; i < Results.size(); i++)
    {
        res_num += Results[i]/Sigma2[i];
        res_den += 1.0/Sigma2[i];
    }
    return res_num/res_den;
}
double VEGAS_Integrator::Get_Error()
{
    double res = 0;
    for (int i = 0; i < Sigma2.size(); i++)
    {
        res += 1.0/Sigma2[i];
    }
    return 1.0/sqrt(res);
}
double VEGAS_Integrator::Get_Chisq()
{
    double Ifinal = Get_Result();
    double chi2 = 0;
    for (int i = 0; i < Results.size(); i++)
    {
        chi2 += pow(Results[i]-Ifinal,2)/Sigma2[i];
    }
    return chi2;
}