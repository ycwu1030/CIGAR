#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;

void VEGAS_Integrator::Set_Verbose(VEGAS_INTEGRATOR_VERBOSE level)
{
    verb = level;
}

void VEGAS_Integrator::Set_Integrand(INTEGRAND integrand, int dim, void* param)
{
    func = integrand;
    N_DIM = dim;
    userdata = param;
    Results.clear();
    Sigma2.clear();
    map = VEGAS_Map(dim);
}

void VEGAS_Integrator::Improve_Grid()
{
    vector<double> yrnd(N_DIM);
    vector<double> y(N_DIM); // Random number between 0 to 1;
    vector<double> x(N_DIM); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; 
    if (verb >= INFO)
    {
        cout<<"==========================================================="<<endl;
        cout<<"| Improving the mapping grid and stratification grid      |"<<endl;
        cout<<"==========================================================="<<endl;
        cout<<"|  Iter  |  N_Eval  |  Sigma [pb]  |  Error [pb]  |  Acc  |"<<endl;
    }
    int iter = 0;
    int NEVAL_START = 10000;
    double dV;
    double Res;
    double Err2;
    int neval;
    int NEVAL_REAL;
    double Jf;
    double Jf2;
    double Ih;
    double Sig2;
    double acc;
    strat.Set_Stratification_System(N_DIM,NEVAL_START);
    while (true)
    {
        // we decide to end the improvement of grid and strata when the accuracy is about 1%
        // Every 5 iteration, we can check the accuracy, and addjust the number of evaluation
        iter++;
        strat.Set_NEVAL(NEVAL_START);
        dV = strat.Get_V_Cubic();
        Res = 0;
        Err2 = 0;
        NEVAL_REAL = 0;
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0;
            Jf2 = 0;
            neval = strat.Get_NH(inc);
            NEVAL_REAL += neval;
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
                map.Accumulate_Weight(y,f_eval);
                strat.Accumulate_Weight(inc,f_eval*Jac);
                Jf += f_eval*Jac;
                Jf2 += pow(f_eval*Jac,2);
            }
            Ih = Jf/neval*dV;
            Sig2 = Jf2/neval*dV*dV - pow(Jf/neval*dV,2);
            Res += Ih;
            Err2 += Sig2/neval;
        }
        map.Update_Map();
        strat.Update_DH();
        acc = sqrt(Err2)/Res;
        if (verb >= INFO)
        {
            cout<<"| "<<iter<<" | "<<NEVAL_REAL<<" | "<<Res<<" | "<<sqrt(Err2)<<" | "<<acc<<" |"<<endl;
        } 
        if (acc < 0.01 && iter >= 5)
        {
            break;
        }
        if (iter%5==0)
        {
            NEVAL_START = NEVAL_START + 5000;//* sqrt(acc/0.01);
        }
    }
    if (verb >= INFO)
    {
        cout<<"==========================================================="<<endl;
    }
}
void VEGAS_Integrator::Integration(double eps_rel, double eps_abs)
{
    // We try to reach either relative error (eps_rel) or absolute error (eps_abs)
    // But we also need to make sure chi2 is not bigger than the iteration numbers
    vector<double> yrnd(N_DIM);
    vector<double> y(N_DIM); // Random number between 0 to 1;
    vector<double> x(N_DIM); // The argument for integrand;
    double f_eval; // evaluated integrand value;
    double Jac; // The Jacobian from y to x;
    int NEVAL_START = 50000;
    double dV = strat.Get_V_Cubic();
    int iter = 0;
    double Res;
    double Err;
    double Chi2;
    int neval;
    int NEVAL_REAL;
    double Jf;
    double Jf2;
    double Ih;
    double Sig2;
    double acc;
    if (verb >= INFO)
    {
        cout<<"==============================================================="<<endl;
        cout<<"| Fixing the mapping grid, still improve strata then Integral |"<<endl;
        cout<<"==============================================================="<<endl;
        cout<<"|  Iter  |  N_Eval  |  Sigma [pb]  |  Error [pb]  |  Acc  |"<<endl;
    }
    while (true)
    {
        iter++;
        strat.Set_NEVAL(NEVAL_START);
        Results.push_back(0);
        Sigma2.push_back(0);
        NEVAL_REAL = 0;
        for (int inc = 0; inc < strat.Get_NHYPERCUBICS(); inc++)
        {
            Jf = 0;
            Jf2 = 0;
            neval = strat.Get_NH(inc);
            NEVAL_REAL += neval;
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
                strat.Accumulate_Weight(inc,f_eval*Jac);
                Jf += f_eval*Jac;
                Jf2 += pow(f_eval*Jac,2);
            }
            Ih = Jf/neval*dV;
            Sig2 = Jf2/neval*dV*dV - pow(Jf/neval*dV,2);
            Results[Results.size()-1] += Ih;
            Sigma2[Sigma2.size()-1] += Sig2/neval;
        }
        strat.Update_DH();
        acc = sqrt(Sigma2[Sigma2.size()-1])/Results[Results.size()-1];
        if (verb >= INFO)
        {
            cout<<"| "<<iter<<" | "<<NEVAL_REAL<<" | "<<Results[Results.size()-1]<<" | "<<sqrt(Sigma2[Sigma2.size()-1])<<" | "<<acc<<" |"<<endl;
        }
        if (iter%5==0)
        {
            // Every 5 iteration, we check whether we fullfil the condition
            Res = Get_Result();
            Err = Get_Error();
            Chi2 = Get_Chisq();
            acc = Err/Res;
            if (verb >= INFO)
            {
                cout<<"| Summary of Last 5 Iter: | Res = "<< Res <<" | Err = "<< Err <<" | Chi2 = "<<Chi2<<" |"<<endl;
            }
            if ( (acc < eps_rel || Err < eps_abs) && Chi2/5.0 < 1.0 )
            {
                break;
            }
            if (Chi2/5.0 < 1.0)
            {
                NEVAL_START = NEVAL_START + 5000;//* sqrt(acc/eps_rel);
                Results.clear();
                Sigma2.clear();
                continue;
            }
            if (Chi2/5.0 > 1.0)
            {
                NEVAL_START += 5000;
                Results.clear();
                Sigma2.clear();
                continue;
            }
        }
    }
    if (verb >= INFO)
    {
        cout<<"==========================================================="<<endl;
        cout<<"Summary: "<<endl;
        cout<<"Result: "<<Get_Result()<<"  Error: "<<Get_Error()<<"  Chi2: "<<Get_Chisq()<<endl;
        cout<<"==========================================================="<<endl;
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