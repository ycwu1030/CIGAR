#include "VEGAS_Stratify.h"
#include <cmath>
#include <iostream>

using namespace std;

void VEGAS_Stratify::Set_Dimension(int ndim)
{
    N_DIM = ndim;
    Reset_Storage();
}
void VEGAS_Stratify::Set_NEVAL(int NEVAL_EXP)
{
    N_EVALUATES_EXPECTED = NEVAL_EXP;
}
// void VEGAS_Stratify::Set_Stratification_System(int ndim, int NEVAL_TRAIN)
// {
//     N_DIM = ndim;
//     N_EVALUATES_TRAINED = NEVAL_TRAIN;
//     Reset_Storage();
// }
void VEGAS_Stratify::Reset_Storage()
{
    // N_STRAT = floor(pow(N_EVALUATES_TRAINED/4.0,1.0/N_DIM));
    
    N_HYPERCUBICS = pow(N_STRAT,N_DIM);
    if (N_HYPERCUBICS > N_HYPERCUBICS_MAX || N_DIM > 9) // if N_DIM too large, N_HYPERCUBICS will exceed the MAXIMUM number an integer can store
    {
        N_STRAT = floor(pow(N_HYPERCUBICS_MAX,1.0/N_DIM));
        N_HYPERCUBICS = pow(N_STRAT,N_DIM);
    }
    
    
    V_cubic = pow(1.0/N_STRAT, N_DIM);
    JF2 = vector<double>(N_HYPERCUBICS,0);
    JF  = vector<double>(N_HYPERCUBICS,0);
    Counts = vector<double>(N_HYPERCUBICS,0);
    dh  = vector<double>(N_HYPERCUBICS,1.0/N_HYPERCUBICS);
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
void VEGAS_Stratify::Accumulate_Weight(int index, double weight)
{
    // This weight is J*f;
    JF2[index] += weight*weight;
    JF[index] += weight;
    Counts[index] += 1;
}
void VEGAS_Stratify::Update_DH()
{
    double d_sum = 0;
    double d_tmp;
    for (int i = 0; i < N_HYPERCUBICS; i++)
    {
        d_tmp = V_cubic*V_cubic/Counts[i]*JF2[i] - pow(V_cubic/Counts[i]*JF[i],2);
        dh[i] = pow(d_tmp,beta);
        d_sum += dh[i];
    }
    for (int i = 0; i < N_HYPERCUBICS; i++)
    {
        dh[i] = dh[i]/d_sum;
    }
}
int VEGAS_Stratify::Get_NH(int index)
{
    int nh = dh[index]*N_EVALUATES_EXPECTED;
    return nh<2?2:nh;
}
