#ifndef VEGAS_MAP_H
#define VEGAS_MAP_H
#include <vector>
#include <string>

// * This is the grid map from y to x:
// * y: the variable upon which we generate uniformly distributed random numbers
// * x: the integral variable, the upper and lower limit set to be 0 and 1 (The mapping from true integrate limits to [0,1] should be done by users in the integrand)
// * The function of this class:
// * 1. Keep track the mapping between y to x
// * 2. Keep the Jacobian from y to x
// * 3. Take care of the grid map improvements
class VEGAS_Map
{
private:
    int N_DIM;
    int N_INTERVALS;
    int N_EDGES; // N_INTERVALS + 1;
    double alpha; // The parameter control the smooth of weight
    
    std::vector<std::vector<double> > x_edges; // The edges in x, size = N_DIM x N_EDGES;
    std::vector<std::vector<double> > dx_steps; // The step for each interval, size = N_DIM x N_INTERVALS;

    std::vector<std::vector<double> > x_edges_last; // The edges in x, size = N_DIM x N_EDGES;
    std::vector<std::vector<double> > dx_steps_last; // The step for each interval, size = N_DIM x N_INTERVALS;
    
    std::vector<std::vector<double> > weights; // The weight in each interval, used to improve the grid map, size = N_DIM x N_INTERVALS;
    std::vector<std::vector<double> > counts; // Count the numbers of random numbers in specific interval, size = N_DIM x N_INTERVALS;
    
    std::vector<std::vector<double> > smoothed_weights; // Smoothed weights, also renormalized, size = N_DIM x  N_INTERVALS
    std::vector<double> summed_weights; // The all summed smoothed weights, size = N_DIM
    std::vector<double> delta_weights; // The step for weights, size = N_DIM

    void Smooth_Weight();
    void Reset_Weight();
    void Checking_Weight();
    std::vector<double> average_weight; // size = N_DIM
    std::vector<double> std_weight; // size = N_DIM

public:
    VEGAS_Map();
    VEGAS_Map(int NDIM);
    VEGAS_Map(int NDIM, int Intervals);
    ~VEGAS_Map(){};

    void Reset_Map();
    void Set_alpha(double alp){alpha = alp;};
    void Accumulate_Weight(std::vector<double> y, double f); // f is the integrand, no other manupulation
    void Update_Map();

    int Get_N_Interval(){return N_INTERVALS;}
    std::vector<int> Get_Interval_ID(std::vector<double> y);
    std::vector<double> Get_Interval_Offset(std::vector<double> y);

    std::vector<double> Get_X(std::vector<double> y);
    double Get_Jac(std::vector<double> y);

    // void Dump_Edges(std::string filename);
    void Print_Edges();
    void Print_Weights();

    double Checking_Map();
};



#endif // VEGAS_MAP_H