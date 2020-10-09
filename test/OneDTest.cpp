#include <iostream>

using namespace std;

extern double VEGAS1D(double xmin, double xmax, double &error);

int main(int argc, char const *argv[])
{
    double result;
    double err;
    for (double cut = 0.001; cut <= 0.11; cut+=0.01)
    {
        result = VEGAS1D(-1+cut,1-cut,err);
        cout<<cut<<"\t"<<result<<endl;
    }
    
    return 0;
}
