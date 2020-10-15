#include <iostream>
#include <cmath>

using namespace std;

extern double VEGAS1D(double sqrts, double &error);

int main(int argc, char const *argv[])
{
    double result;
    double err;
    double energies[37] = {1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000,22000,24000,26000,28000,30000};
    for (int i = 0; i < 37; i++)
    {
        result = VEGAS1D(energies[i],err);
        cout<<energies[i]<<"\t"<<result/pow(energies[i],2)<<endl;
    }
    
    return 0;
}
