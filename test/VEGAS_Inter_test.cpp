#include "VEGAS_Integrator.h"
#include <iostream>

using namespace std;

double func_weight(vector<double> x, void* param)
{
    double dx = *((double *)param);
    double xmin = -1.0 + dx;
    double xmax = 1.0 - dx;
    double x_true = x[0]*2.0-1.0;
    if (x_true < xmin || x_true > xmax)
    {
        return 0;
    }
    return (1.0+x_true*x_true)/(1.0-x_true*x_true)*2;
}
double func_WW2hh(vector<double> x, void *param)
{
    double *par = (double *)param;
    double s = par[0]*par[0];
    double MH = 125.0;
    double MW = 80.385;
    double MH2 = MH*MH;
    double MW2 = MW*MW;
    double cth = x[0]*2.0-1.0;
    double cth2 = cth*cth;
    double num = 8*pow(MH2,3)*(cth2*s-2*MW2+s)+2*MH2*MH2*(16*(2*cth2+1)*MW2*MW2-40*cth2*MW2*s+(cth2-3)*s*s) + MH2*s*(16*(cth2-3)*MW2*MW2 + 4*(7*cth2+4)*MW2*s-(cth2-1)*s*s) - 8*(cth2-2)*MW2*MW2*s*s - 2*(cth2+3)*MW2*s*s*s;
    double den = 16*MW2*MW2*pow(s-MH2,2)*pow(-4*MH2*(cth2*(4*MW2-s)+s)+s*(cth2*(4*MW2-s)+s)+4*MH2*MH2,2);
    return pow(num,2)/den;
}
double func_4D(vector<double> x, void *param)
{
    double r1[4] = {0.33,0.5,0.5,0.5};
    double r2[4] = {0.67,0.5,0.5,0.5};
    double x1=0;
    double x2=0;
    for (int i = 0; i < 4; i++)
    {
        x1 += -100*pow(x[i]-r1[i],2);
        x2 += -100*pow(x[i]-r2[i],2);
    }
    return exp(x1) + exp(x2);
}
int main(int argc, char const *argv[])
{
    VEGAS_Integrator inter;
    
    for (double dx = 0.02; dx < 0.31; dx += 0.02)
    {
        inter.Set_Integrand(func_weight, 1, &dx);
        inter.Improve_Grid();
        inter.Integration();
        cout<<"dx: "<<dx<<" res: "<<inter.Get_Result()<<" err: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<endl;
    }
    double energies[37] = {1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000,22000,24000,26000,28000,30000};
    for (int i = 0; i < 37; i++)
    {
        inter.Set_Integrand(func_WW2hh, 1, &energies[i]);
        inter.Improve_Grid();
        inter.Integration();
        cout<<"Ecm: "<<energies[i]<<" res: "<<inter.Get_Result()/pow(energies[i],2)<<" err: "<<inter.Get_Error()/pow(energies[i],2)<<" chi2: "<<inter.Get_Chisq()<<endl;
    }
    inter.Set_Integrand(func_4D,4,NULL);
    inter.Improve_Grid();
    inter.Integration();
    cout<<"Result: "<<inter.Get_Result()<<" Error: "<<inter.Get_Error()<<" chi2: "<<inter.Get_Chisq()<<endl;
    return 0;
}
