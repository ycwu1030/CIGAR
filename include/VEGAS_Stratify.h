#ifndef VEGAS_STRATIFY_H
#define VEGAS_STRATIFY_H

#include <vector>

class VEGAS_Stratify
{
private:
    int N_DIM;
    int N_STRAT;
    double beta;
    double V_cubic;
    std::vector<double> JF2; // size = (N_STRAT)^(N_DIM)
    std::vector<double> JF; // size = (N_STRAT)^(N_DIM)
    std::vector<int> nh; // size = (N_STRAT)^(N_DIM)
    std::vector<double> dh; // size = (N_STRAT)^(N_DIM)
    int N_EVALUATES_EXPECTED;
    int N_EVALUATES_ACTUALLY;
    int N_HYPERCUBICS;

    void Reset_Storage();
    std::vector<int> Get_Indices(int index);

public:
    VEGAS_Stratify(){N_DIM = 1; N_EVALUATES_EXPECTED = 10000; beta = 0.75;};
    ~VEGAS_Stratify(){};

    void Set_Stratification_System(int N_DIM, int NEVAL);
    void Set_Dimension(int N_DIM);
    void Set_NEVAL(int NEVAL);
    void Accumulate_Weights(int index, double weight);
    void Update_NH();
    std::vector<double> Get_Y(int index, std::vector<double> random_uni);
    int Get_NH(int index){}; // Get the expected number of events in each hypercubic.
};


#endif //VEGAS_STRATIFY_H