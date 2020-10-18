#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;


void VEGAS_Integrator::Set_Integrand(INTEGRAND integrand, int dim, void* param)
{
    func = integrand;
    N_DIM = dim;
    userdata = param;
    Results.clear();
    Sigma2.clear();
    map = VEGAS_Map(dim);
}

void VEGAS_Integrator::Improve_Grid(int Iter, int Neval)
{
    vector<double> yrnd(N_DIM);
    vector<double> y(N_DIM); // Random number between 0 to 1;
    vector<double> x(N_DIM); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; 
    strat.Set_Stratification_System(N_DIM,Neval);
    strat.Set_NEVAL(Neval);
    for (int it = 0; it < Iter; it++)
    {
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            for (int ne = 0; ne < strat.Get_NH(inc); ne++)
            {
                for (int i_dim = 0; i_dim < N_DIM; i_dim++)
                {
                    yrnd[i_dim] = dist(rng);
                }
                y = strat.Get_Y(inc,yrnd);
                x = map.Get_X(y);
                f_eval = func(x,userdata);
                Jac = map.Get_Jac(y);
                map.Accumulate_Weight(y,f_eval);
                strat.Accumulate_Weight(inc,f_eval*Jac);
            }
        }
        map.Update_Map();
        strat.Update_DH();
    }
}
void VEGAS_Integrator::Integration(int Iter, int Neval)
{
    vector<double> yrnd(N_DIM);
    vector<double> y(N_DIM); // Random number between 0 to 1;
    vector<double> x(N_DIM); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; // The Jacobian from y to x;
    strat.Set_NEVAL(Neval);
    double dV = strat.Get_V_Cubic();
    for (int it = 0; it < Iter; it++)
    {
        Results.push_back(0);
        Sigma2.push_back(0);
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            double Jf = 0;
            double Jf2 = 0;
            int neval = strat.Get_NH(inc);
            for (int ne = 0; ne < neval; ne++)
            {
                for (int i_dim = 0; i_dim < N_DIM; i_dim++)
                {
                    yrnd[i_dim] = dist(rng);
                }
                y = strat.Get_Y(inc,yrnd);
                x = map.Get_X(y);
                f_eval = func(x,userdata);
                Jac = map.Get_Jac(y);
                Jf += f_eval*Jac;
                Jf2 += pow(f_eval*Jac,2);
            }
            double Ih = Jf/neval*dV;
            double Sig2 = Jf2/neval*dV*dV - pow(Jf/neval*dV,2);
            Results[Results.size()-1] += Ih;
            Sigma2[Sigma2.size()-1] += Sig2/neval;
        }
        // Results.push_back(res/Neval);
        // Sigma2.push_back((sig2/Neval-pow(res/Neval,2))/(Neval-1.0));
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