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
    std::vector<double> Counts; // size = (N_STRAT)^(N_DIM)
    std::vector<double> dh; // size = (N_STRAT)^(N_DIM)
    // int N_EVALUATES_TRAINED; // The evaluates number used to train the stratification
    int N_EVALUATES_EXPECTED;
    int N_HYPERCUBICS;
    int N_HYPERCUBICS_MAX;

    void Reset_Storage();
    std::vector<int> Get_Indices(int index);

public:
    VEGAS_Stratify(){N_DIM = 1; N_STRAT = 10; beta = 0.75; N_HYPERCUBICS_MAX = 10000;};
    ~VEGAS_Stratify(){};

    // void Set_Stratification_System(int N_DIM, int NEVAL_TRAIN);
    void Set_Dimension(int N_DIM);
    void Set_NEVAL(int NEVAL_EXP);
    void Accumulate_Weight(int index, double weight);
    void Update_DH();
    std::vector<double> Get_Y(int index, std::vector<double> random_uni);
    int Get_NHYPERCUBICS(){return N_HYPERCUBICS;};
    int Get_NH(int index); // Get the expected number of events in each hypercubic.
    double Get_V_Cubic(){return V_cubic;}
};


#endif //VEGAS_STRATIFY_H