#include "VEGAS_map.h"
#include <random> // The random number generator and distributions. c++11
#include <vector>

typedef double (*INTEGRAND)(double x, void *param);
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
    INTEGRAND func;
    void* userdata;

    VEGAS_Map map;

    std::mt19937 rng; // Mersenne twister random number engine
    URD dist; // uniform distribution in double in [0.0, 1.0)

    std::vector<double> Results;
    std::vector<double> Sigma2;

public:
    VEGAS_Integrator(){};
    ~VEGAS_Integrator(){};

    void Set_Integrand(INTEGRAND integrand, void* param);
    void Improve_Grid(int Iter = 5, int Neval = 10000);
    void Integration(int Iter = 5, int Neval = 50000);
    
    
    double Get_Result();
    double Get_Error();
    double Get_Chisq();

};


