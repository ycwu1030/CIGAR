// Try to follow VEGAS's algorithm to integrate over (1+x^2)/(1-x^2) from very close to -1 to very close to 1
#include <vector>
#include <random>
#include <algorithm>
#include <TH1F.h>
#include <iostream>

using namespace std;
int FindPosition(const vector<double> &xi, const double x)
{
    // Do the binary search, try to reach O(log(N))
    int NSize = xi.size();
    if (x<xi[0]||x>xi[NSize-1]) return -1;
    int icur = NSize/2;
    while (true)
    {
        if (xi[icur] <= x && x < xi[icur+1])
        {
            return icur;
        }
        else if (x<xi[icur])
        {
            icur /= 2;
        }
        else
        {
            icur = (NSize + icur)/2;
        }
    }
}
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
        partition[cur-1] += 1;
    }
    // cout<<"\t\tPartition: ";
    // for (int i = 0; i < N; i++)
    // {
    //     cout<<" "<<partition[i];
    // }
    // cout<<endl;
    vector<double> xori = xi;
    int current_bigN=0;
    int current_littleN=0;
    for (int i = 1; i < N; i++)
    {
        // cout<<"i="<<i<<endl;
        current_littleN += partition[i-1];
        // cout<<"\t\tCurrent BIG: "<<current_bigN<<endl;
        // cout<<"\t\tCurrent lit: "<<current_littleN<<endl;
        // cout<<"\t\tMIS: "<<mis[current_bigN]<<endl;
        while (current_littleN > mis[current_bigN])
        {
            current_littleN-=mis[current_bigN];
            current_bigN++;
            // cout<<"\t\t\tCurrent BIG: "<<current_bigN<<endl;
            // cout<<"\t\t\tCurrent lit: "<<current_littleN<<endl;
            // cout<<"\t\t\tMIS: "<<mis[current_bigN]<<endl;
        }
        xi[i] = xori[current_bigN] + current_littleN*(xori[current_bigN+1]-xori[current_bigN])/mis[current_bigN];
        // cout<<"xi["<<i<<"]="<<xi[i]<<endl;
    }
    return;
}
double Integrand(double x, void *param)
{
    double *par = (double *)param;
    double xmin = par[0];
    double xmax = par[1];
    // double x_true = x*(xmax-xmin)+xmin;
    double x_true = x*2.0-1.0;
    if (x_true < xmin || x_true > xmax)
    {
        return 0;
    }
    return (1.0+x_true*x_true)/(1.0-x_true*x_true)*2;//(xmax-xmin);
}

double VEGAS1D(double xmin, double xmax, double &error)
{
    // Constants declear
    const int N = 80;  // The fixed number of increments from xmin to xmax
    const int K = 1000; // The number used to refined the increments
    const double alpha = 1.5;
    const double EPS_ABS = 1e-9;
    const double EPS_REL = 1e-3;
    int NM = 2000;


    // Initiation
    vector<double> xi_N_Zero2One(N+1);
    vector<double> weight_N(N);
    vector<int> mis(N);
    for (int i = 0; i < N+1; i++)
    {
        xi_N_Zero2One[i] = i*1.0/N;
        if (i>0)
        {
            weight_N[i-1]=1.0/N;
        }
    }

    default_random_engine generator;
    double x_rand;
    double f_integrand;
    double xrange[2] = {xmin,xmax};
    bool first_time = true;
    double result_accumulate=0;
    double result2_accumulate=0;
    double result_cur=0;
    double err_cur=0;
    double err_last=0;
    double f2_results_cur;
    double abs_f_result;
    double px;
    int iter = 0;
    while (true)
    {
        iter++;
        // Start the interation
        piecewise_constant_distribution<double> dist(xi_N_Zero2One.begin(),xi_N_Zero2One.end(),weight_N.begin());
        // vector<double> h1(N,0);
        // vector<double> h2(N,0);
        // vector<double> h3(N,0);
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
        err_last = err_cur;
        result_cur = 0;
        f2_results_cur = 0;
        abs_f_result = 0;
        for (int i = 0; i < N; i++)
        {
            px = 1.0/(xi_N_Zero2One[i+1]-xi_N_Zero2One[i])/N;
            result_cur += h2->GetBinContent(i+1)/px;
            f2_results_cur += h3->GetBinContent(i+1)/px/px;
            abs_f_result += h1->GetBinContent(i+1)*(xi_N_Zero2One[i+1]-xi_N_Zero2One[i]);
        }
        result_cur /= NM;
        f2_results_cur /= NM;
        err_cur = sqrt((f2_results_cur - pow(result_cur,2))/(NM-1.0));

        result_accumulate += result_cur*result_cur/err_cur/err_cur;
        result2_accumulate += pow(result_cur,3)/err_cur/err_cur;
        result_cur = result2_accumulate/result_accumulate;
        err_cur = result_cur/sqrt(result_accumulate);

        // cout<<"\tIter: "<<iter<<" Result: "<<result_cur<<" err: "<<err_cur<<" errlast: "<<err_last<<" NM: "<<NM<<endl;

        if (iter > 1)
        {
            if (err_cur/err_last > 0.95)
            {
                NM *= 4;
                NM = NM>200000?200000:NM;
            }
        }
        

        
        if (err_cur < EPS_ABS || err_cur/result_cur < EPS_REL)
        {
            delete h1;
            delete h2;
            delete h3;
            break;
        }
        
        // Improve the grid
        int MTOTAL = 0;
        // cout<<"\t\tINCREMENTS: ";
        for (int i = 0; i < N; i++)
        {
            double portion = h1->GetBinContent(i+1)*(xi_N_Zero2One[i+1]-xi_N_Zero2One[i])/abs_f_result;
            portion=portion==0?1e-15:portion;
            mis[i] = floor(K*pow((portion-1.0)/log(portion),alpha)+1);
            MTOTAL += mis[i];
            // cout<<" "<<mis[i];
        }
        // cout<<endl;
        ImproveGrid(N,mis,MTOTAL,xi_N_Zero2One);
        // cout<<"\t\tGRID: ";
        // for (int i = 0; i < N+1; i++)
        // {
        //     cout<<" "<<xi_N_Zero2One[i];
        // }
        // cout<<endl;
        for (int i = 1; i < N+1; i++)
        {
            weight_N[i-1]=1.0/N;
        }
        delete h1;
        delete h2;
        delete h3;
    }
    
    error = err_cur;
    return result_cur;
    

}

// int main(int argc, char const *argv[])
// {
//     double result;
//     double err;
//     for (double cut = 0.01; cut <= 0.21; cut+=0.01)
//     {
//         result = VEGAS1D(-1+cut,1-cut,err);
//         cout<<cut<<"\t"<<result<<endl;
//     }
    
//     return 0;
// }