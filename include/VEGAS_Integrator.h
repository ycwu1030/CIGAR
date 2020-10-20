#ifndef VEGAS_INTEGRATOR_H
#define VEGAS_INTEGRATOR_H

#include "VEGAS_map.h"
#include "VEGAS_Stratify.h"
#include <random> // The random number generator and distributions. c++11
#include <vector>

typedef double (*INTEGRAND)(std::vector<double> x, void *param);
using URD=std::uniform_real_distribution<double>;

enum VEGAS_INTEGRATOR_VERBOSE
{
    NONE = 0,
    INFO = 1,
    DEBUG = 2,
    ALL = 3
};

class VEGAS_Integrator
{
private:
    VEGAS_INTEGRATOR_VERBOSE verb;

    INTEGRAND func;
    int N_DIM;
    void* userdata;

    VEGAS_Map map;
    VEGAS_Stratify strat;

    std::mt19937 rng; // Mersenne twister random number engine
    URD dist; // uniform distribution in double in [0.0, 1.0)

    std::vector<double> Results;
    std::vector<double> Sigma2;


public:
    VEGAS_Integrator(){verb = INFO;};
    ~VEGAS_Integrator(){};

    void Set_Verbose(VEGAS_INTEGRATOR_VERBOSE level);

    void Set_Integrand(INTEGRAND integrand, int dim, void* param);
    void Improve_Grid();
    void Integration(double eps_rel = 1e-3, double eps_abs = 1e-9);
    
    
    double Get_Result();
    double Get_Error();
    double Get_Chisq();

};


#endif //VEGAS_INTEGRATOR_H