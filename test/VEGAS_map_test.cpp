#include "VEGAS_map.h"
#include <random>
#include <iostream>

using namespace std;

double func_weight(vector<double> x)
{
    double dx = 0.23;
    double xmin = -1.0 + dx;
    double xmax = 1.0 - dx;
    double x_true = x[0]*2.0-1.0;
    if (x_true < xmin || x_true > xmax)
    {
        return 0;
    }
    return (1.0+x_true*x_true)/(1.0-x_true*x_true)*2;
}
// double func_weight(double x)
// {
//     double smin = 200*200;
//     double smax = 1000*1000;
//     double mm = 800;
//     double gamma = 10;
//     double s = x*(smax-smin)+smin;
//     return 1.0/(pow(s-mm*mm,2)+mm*mm*gamma*gamma);
// }
int main(int argc, char const *argv[])
{
    VEGAS_Map GM;
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0,1);
    for (int i = 0; i < 50; i++)
    {
        cout<<"=========>"<<endl;
        // GM.Print_Edges();
        vector<double> y(1);
        vector<double> x(1);
        for (int j = 0; j < 10000; j++)
        {
            y[0] = distribution(generator);
            x = GM.Get_X(y);
            double weight = func_weight(x);
            GM.Accumulate_Weight(y,weight);
        }
        // GM.Print_Weights();
        GM.Update_Map();
        GM.Print_Edges();
        cout<<"<========="<<endl;
    }
    
    return 0;
}
