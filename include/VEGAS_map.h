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
    int N_INTERVALS;
    int N_EDGES; // N_INTERVALS + 1;
    double alpha; // The parameter control the smooth of weight
    
    std::vector<double> x_edges; // The edges in x, size = N_EDGES;
    std::vector<double> dx_steps; // The step for each interval, size = N_INTERVALS;
    
    std::vector<double> weights; // The weight in each interval, used to improve the grid map, size = N_INTERVALS;
    std::vector<double> counts; // Count the numbers of random numbers in specific interval
    
    std::vector<double> smoothed_weights; // Smoothed weights, also renormalized.
    double summed_weights; // The all summed smoothed weights
    double delta_weights; // The step for weights

    void Smooth_Weight();
    void Reset_Weight();
    void Checking_Weight();
    double average_weight;
    double std_weight;

public:
    VEGAS_Map();
    VEGAS_Map(int Intervals);
    ~VEGAS_Map(){};

    void Reset_Map();
    void Accumulate_Weight(double y, double f); // f is the integrand, no other manupulation
    void Update_Map();

    int Get_N_Interval(){return N_INTERVALS;}
    int Get_Interval_ID(double y);
    double Get_Interval_Offset(double y);

    double Get_X(double y);
    double Get_Jac(double y);

    // void Dump_Edges(std::string filename);
    void Print_Edges();
    void Print_Weights();
};



#endif // VEGAS_MAP_H