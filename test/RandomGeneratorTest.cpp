#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TH1F.h>
#include <TCanvas.h>

using namespace std;


int main(int argc, char const *argv[])
{
    const int N = 10;
    vector<double> xi;
    vector<double> wi(N);
    xi.push_back(0.0);
    xi.push_back(0.06);
    xi.push_back(0.39);
    xi.push_back(0.43);
    xi.push_back(0.52);
    xi.push_back(0.57);
    xi.push_back(0.62);
    xi.push_back(0.81);
    xi.push_back(0.88);
    xi.push_back(0.93);
    xi.push_back(1.0);

    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    double xmintmp = 0.0;
    double xmaxtmp = 1.0;
    // for (int i = 1; i <= N-1; i++)
    // {
    //     double tmp = distribution(generator);
    //     double num = tmp*(xmaxtmp-xmintmp)+xmintmp;
    //     if (abs(xmaxtmp-num) > abs(num-xmintmp))
    //     {
    //         xmintmp = num;
    //     }
    //     else
    //     {
    //         xmaxtmp = num;
    //     }
    //     xi.push_back(num);
    // }
    
    sort(xi.begin(),xi.end());
    for (int i = 0; i < N; i++)
    {
        wi[i] = 1.0/(xi[i+1]-xi[i])/N;
        // cout<<xi[i]<<endl;
    }
    
    const int NM = 100000;
    piecewise_constant_distribution<double> dist(xi.begin(),xi.end(),wi.begin());
    TH1F *h1 = new TH1F("h1","",N,xi.data());
    for (int i = 0; i < NM; i++)
    {
        double num = dist(generator);
        h1->Fill(num);
    }
    TCanvas *c1 = new TCanvas("c1","",800,600);
    h1->Draw();
    c1->SaveAs("./auxillary/random_test.png");

    return 0;
}
