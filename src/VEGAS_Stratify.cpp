#include "VEGAS_Stratify.h"

using namespace std;

void VEGAS_Stratify::Set_Dimension(int ndim)
{
    N_DIM = ndim;
}
void VEGAS_Stratify::Set_NEVAL(int NEVAL)
{
    N_EVALUATES_EXPECTED = NEVAL;
}
void VEGAS_Stratify::Set_Stratification_System(int ndim, int NEVAL)
{
    N_DIM = ndim;
    N_EVALUATES_EXPECTED = NEVAL;
    Reset_Storage();
}
void VEGAS_Stratify::Reset_Storage()
{
    N_STRAT = floor(N_EVALUATES_EXPECTED/4.0,1.0/N_DIM);
    
    N_HYPERCUBICS = pow(N_STRAT,N_DIM);

    V_cubic = pow(1.0/N_STRAT, N_DIM);
    JF2 = vector<double>(N_HYPERCUBICS,0);
    JF  = vector<double>(N_HYPERCUBICS,0);
    nh  = vector<int>(N_HYPERCUBICS,N_EVALUATES_EXPACTED/N_HYPERCUBICS);
    dh  = vector<double>(N_HYPERCUBICS,0);
    
    N_EVALUATES_ACTUALLY = 0;

}
vector<int> VEGAS_Stratify::Get_Indices(int index)
{
    vector<int> res(N_DIM,0);
    int Quotient;
    int tmp = index;
    int Remainder;
    for (int i = 0; i < N_DIM; i++)
    {
        Quotient = tmp/N_STRAT;
        Remainder = tmp - Quotient*N_STRAT;
        res[i] = Remainder;
        tmp = Quotient;
    }
    return res;
}
vector<double> VEGAS_Stratify::Get_Y(int index, vector<double> random_uni)
{
    double dy = 1.0/N_STRAT;
    vector<double> res(N_DIM,0);
    vector<int> ID = Get_Indices(index);
    for (int i = 0; i < N_DIM; i++)
    {
        res[i] = random_uni[i]*dy + ID[i]*dy;
    }
    return res;
}
void VEGAS_Stratify::Accumulate_Weights(int index, double weight)
{
    // This weight is J*f;
    JF2[index] += weight*weight;
    JF[index] += weight;
}
void VEGAS_Stratify::Update_NH()
{
    double d_sum = 0;
    double d_tmp;
    for (int i = 0; i < N_HYPERCUBICS; i++)
    {
        d_tmp = V_cubic*V_cubic/nh[i]*JF2[i] - pow(V_cubic/nh[i]*JF[i],2);
        dh[i] = pow(d_tmp,beta);
        d_sum += dh[i];
    }
    int nh_tmp;
    for (int i = 0; i < N_HYPERCUBICS; i++)
    {
        nh_tmp = N_EVALUATES_EXPECTED * dh[i]/d_sum;
        nh[i] = nh_tmp >= 2?nh_tmp:2;
    }
}
