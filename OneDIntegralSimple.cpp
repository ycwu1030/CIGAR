// Try to follow VEGAS's algorithm to integrate over (1+x^2)/(1-x^2) from very close to -1 to very close to 1
#include <vector>
#include <random>
#include <algorithm>
#include <TH1F.h>
#include <iostream>

using namespace std;
void ImproveGrid(const int N, const vector<int> mis, const int MTOTAL, vector<double> &xi)
{
    int m_each = MTOTAL/N;
    int m_left = MTOTAL - m_each*N;
    vector<int> partition(N,m_each);
    default_random_engine generator;
    uniform_int_distribution<int> distribution (1,N);
    for (int i = 0; i < m_left; i++)
    {
        int cur = distribution(generator);
        partition[cur] += 1;
    }
    for (int i = 0; i < N; i++)
    {
        cout<<partition[i]<<" ";
    }
    cout<<endl;
    vector<double> xori = xi;
    int current_bigN=0;
    int current_littleN=0;
    for (int i = 1; i < N; i++)
    {
        current_littleN += partition[i-1];
        while (current_littleN > mis[current_bigN])
        {
            current_littleN-=mis[current_bigN];
            current_bigN++;
        }
        xi[i] = xori[current_bigN] + current_littleN*(xori[current_bigN+1]-xori[current_bigN])/mis[current_bigN];
    }
}
double Integrand(double x, void *param)
{
    double *par = (double *)param;
    double xmin = par[0];
    double xmax = par[1];
    double x_true = x*(xmax-xmin)+xmin;
    return (1.0+x_true*x_true)/(1.0-x_true*x_true);
}

double VEGAS1D(double xmin, double xmax, double &error)
{
    // Constants declear
    const int N = 80;  // The fixed number of increments from xmin to xmax
    const int K = 1000; // The number used to refined the increments
    const double alpha = 1.5;
    const double EPS_ABS = 1e-9;
    const double EPS_REL = 1e-3;
    const int NM = 2000;


    // Initiation
    vector<double> xi_N_Zero2One(N+1);
    vector<double> weight_N(N);
    for (int i = 0; i < N+1; i++)
    {
        xi_N_Zero2One[i] = i*1.0/N;
        if (i>0)
        {
            weight_N[i-1]=1.0/(xi_N_Zero2One[i]-xi_N_Zero2One[i-1])/N;
        }
    }

    default_random_engine generator;
    double x_rand;
    double f_integrand;
    double xrange[2] = {xmin,xmax};
    bool first_time = true;
    double result_last=0,result_cur=0;
    double err_last=0,err_cur=0;
    double f2_results_cur;
    double px;
    while (true)
    {
        // Start the interation
        piecewise_constant_distribution<double> dist(xi_N_Zero2One.begin(),xi_N_Zero2One.end(),weight_N.begin());
        TH1F *h1 = new TH1F("h1","",N,xi_N_Zero2One.data());
        TH1F *h2 = new TH1F("h2","",N,xi_N_Zero2One.data());
        TH1F *h3 = new TH1F("h3","",N,xi_N_Zero2One.data());
        // Generate the random numbers according the distribution
        for (int i = 0; i < NM; i++)
        {
            x_rand = dist(generator);
            f_integrand = Integrand(x_rand,(void*)xrange);
            h1->Fill(x_rand,abs(f_integrand));
            h2->Fill(x_rand,f_integrand);
            h3->Fill(x_rand,pow(f_integrand,2));
        }

        // Estimate the integral
        result_last = result_cur;
        err_last = err_cur;
        result_cur = 0;
        f2_results_cur = 0;
        for (int i = 1; i <= N; i++)
        {
            px = 1.0/(xi_N_Zero2One[i]-xi_N_Zero2One[i-1])/N;
            result_cur += h2->GetBinContent(i)/px;
            f2_results_cur += h3->GetBinContent(i)/px/px;
        }
        result_cur /= NM;
        f2_results_cur /= NM;
        err_cur = sqrt((f2_results_cur - pow(result_cur,2))/(NM-1.0));
        err_cur = max(err_cur,result_cur-result_last);

        if (err_cur < EPS_ABS || err_cur/result_cur < EPS_REL)
        {
            break;
        }
        
        // Improve the grid
        
        

    }
    
    return 0;
    

}
