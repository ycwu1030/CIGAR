#include <vector>
#include <random>
#include <iostream>

using namespace std;

extern void ImproveGrid(const int N, const vector<int> mis, const int MTOTAL, vector<double> &xi);

int main(int argc, char const *argv[])
{
    const int N = 10;
    
    vector<double> xi = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    vector<int> mis;
    int MTOTAL = 0;
    default_random_engine generator;
    uniform_int_distribution<int> distribution(1,20);
    for (int i = 0; i < N; i++)
    {
        int mcur = distribution(generator);
        MTOTAL += mcur;
        mis.push_back(mcur);
    }
    ImproveGrid(N,mis,MTOTAL,xi);

    for (int i = 0; i < N; i++)
    {
        cout<<mis[i]<<"\t"<<xi[i]<<endl;
    }
    cout<<"  \t"<<xi[N]<<endl;
    return 0;
}
