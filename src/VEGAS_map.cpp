#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>

using namespace std;

VEGAS_Map::VEGAS_Map()
{
    N_INTERVALS = 1000;
    N_EDGES = N_INTERVALS + 1;
    alpha = 1.5;
    Reset_Map();
}
VEGAS_Map::VEGAS_Map(int Intervals)
{
    N_INTERVALS = Intervals;
    N_EDGES = N_INTERVALS + 1;
    alpha = 1.5;
    Reset_Map();
}
void VEGAS_Map::Reset_Map()
{
    x_edges.clear();
    dx_steps.clear();
    double step_tmp = 1.0/N_INTERVALS;
    for (int i = 0; i < N_EDGES; i++)
    {
        x_edges.push_back(i*step_tmp);
        if (i>0)
        {
            dx_steps.push_back(x_edges[i]-x_edges[i-1]);
        }
    }
    weights = vector<double>(N_INTERVALS,0);
    counts = vector<double>(N_INTERVALS,0);
    smoothed_weights = vector<double>(N_INTERVALS,0);
}
void VEGAS_Map::Reset_Weight()
{
    weights = vector<double>(N_INTERVALS,0);
    counts = vector<double>(N_INTERVALS,0);
}

int VEGAS_Map::Get_Interval_ID(double y)
{
    return floor(y*N_INTERVALS);
}
double VEGAS_Map::Get_Interval_Offset(double y)
{
    return y*N_INTERVALS - Get_Interval_ID(y);
}
double VEGAS_Map::Get_X(double y)
{
    int id = Get_Interval_ID(y);
    double offset = Get_Interval_Offset(y);
    return x_edges[id] + dx_steps[id]*offset;
}
double VEGAS_Map::Get_Jac(double y)
{
    int id = Get_Interval_ID(y);
    return N_INTERVALS*dx_steps[id];
}
void VEGAS_Map::Accumulate_Weight(double y, double f)
{
    // f is the value of integrand!
    int id = Get_Interval_ID(y);
    weights[id] += pow(f*Get_Jac(y),2);
    counts[id] += 1;
}
void VEGAS_Map::Smooth_Weight()
{
    for (int i = 0; i < weights.size(); i++)
    {
        if (counts[i]!=0)
        {
            weights[i]/=counts[i];
        }
    }
    double d_tmp;
    double d_sum = accumulate(weights.begin(),weights.end(),0.0);
    summed_weights = 0;
    for (int i = 0; i < N_INTERVALS; i++)
    {
        if (i==0)
        {
            d_tmp = (7.0*weights[0]+weights[1])/(8.0*d_sum);
            if (d_tmp == 0)
            {
                d_tmp = 0;
            }
            else
            {
                d_tmp = pow((d_tmp-1.0)/log(d_tmp),alpha);
            }
        }
        else if (i==N_INTERVALS-1)
        {
            d_tmp = (weights[N_INTERVALS-2]+7.0*weights[N_INTERVALS-1])/(8.0*d_sum);
            if (d_tmp == 0)
            {
                d_tmp = 0;
            }
            else
            {
                d_tmp = pow((d_tmp-1.0)/log(d_tmp),alpha);
            }
        }
        else
        {
            d_tmp = (weights[i-1] + 6.0*weights[i] + weights[i+1])/(8.0*d_sum);
            if (d_tmp == 0)
            {
                d_tmp = 0;
            }
            else
            {
                d_tmp = pow((d_tmp-1.0)/log(d_tmp),alpha);
            }
        }
        smoothed_weights[i]=d_tmp;
        summed_weights += d_tmp;
    }
    delta_weights = summed_weights/N_INTERVALS;
}
void VEGAS_Map::Update_Map()
{
    Smooth_Weight();
    vector<double> x_edges_old = x_edges;
    vector<double> dx_steps_old = dx_steps;
    int current_old_interval = 0;
    int current_new_interval = 1;
    double d_accu = 0;
    while (true)
    {
        d_accu += delta_weights;
        while (d_accu > smoothed_weights[current_old_interval])
        {
            d_accu -= smoothed_weights[current_old_interval];
            current_old_interval++;
        }
        x_edges[current_new_interval] = x_edges_old[current_old_interval] + d_accu/smoothed_weights[current_old_interval]*dx_steps_old[current_old_interval];
        dx_steps[current_new_interval-1] = x_edges[current_new_interval] - x_edges[current_new_interval-1];
        current_new_interval++;
        if (current_new_interval >= N_INTERVALS)
        {
            break;
        }
        
    }
    dx_steps[N_INTERVALS-1] = x_edges[N_EDGES-1] - x_edges[N_EDGES-2];
    Reset_Weight();
}
void VEGAS_Map::Checking_Weight()
{
    average_weight = 0;
    for (int i = 0; i < weights.size(); i++)
    {
        average_weight += weights[i];
    }
    average_weight /= weights.size();
    for (int i = 0; i < weights.size(); i++)
    {
        std_weight += pow(weights[i]-average_weight,2);
    }
    std_weight = sqrt(std_weight);// /average_weight;
}
void VEGAS_Map::Print_Edges()
{
    cout<<"Grid Map:";
    for (int i = 0; i < N_EDGES; i++)
    {
        cout<<"\t"<<x_edges[i];
    }
    cout<<endl;
    cout<<"Grid Gap:";
    for (int i = 0; i < N_INTERVALS; i++)
    {
        cout<<"\t"<<dx_steps[i];
    }
    cout<<endl;
}
void VEGAS_Map::Print_Weights()
{
    cout<<"Weights:";
    for (int i = 0; i < N_INTERVALS; i++)
    {
        cout<<"\t"<<weights[i];
    }
    cout<<endl;
    Checking_Weight();
    cout<<"Average: "<<average_weight<<"  STD: "<<std_weight<<" std/ave: "<<std_weight/average_weight/N_INTERVALS<<endl;
}