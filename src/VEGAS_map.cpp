#include "VEGAS_map.h"
#include <cmath>
#include <numeric>
#include <iostream>

using namespace std;

VEGAS_Map::VEGAS_Map()
{
    N_DIM = 1;
    N_INTERVALS = 1000;
    N_EDGES = N_INTERVALS + 1;
    alpha = 0.5;
    Reset_Map();
}
VEGAS_Map::VEGAS_Map(int NDIM)
{
    N_DIM = NDIM;
    N_INTERVALS = 1000;
    N_EDGES = N_INTERVALS + 1;
    alpha = 0.5;
    Reset_Map();
}
VEGAS_Map::VEGAS_Map(int NDIM, int Intervals)
{
    N_DIM = NDIM;
    N_INTERVALS = Intervals;
    N_EDGES = N_INTERVALS + 1;
    alpha = 0.5;
    Reset_Map();
}
void VEGAS_Map::Reset_Map()
{
    x_edges.clear();
    dx_steps.clear();
    double step_tmp = 1.0/N_INTERVALS;
    vector<double> x_edges_tmp;
    vector<double> dx_steps_tmp;
    for (int i = 0; i < N_EDGES; i++)
    {
        x_edges_tmp.push_back(i*step_tmp);
        if (i>0)
        {
            dx_steps_tmp.push_back(x_edges_tmp[i]-x_edges_tmp[i-1]);
        }
    }
    x_edges = vector<vector<double> >(N_DIM,x_edges_tmp);
    dx_steps = vector<vector<double> >(N_DIM,dx_steps_tmp);
    x_edges_last = vector<vector<double> >(N_DIM,x_edges_tmp);
    dx_steps_last = vector<vector<double> >(N_DIM,dx_steps_tmp);
    weights = vector<vector<double> >(N_DIM,vector<double>(N_INTERVALS,0));
    counts = vector<vector<double> >(N_DIM,vector<double>(N_INTERVALS,0));
    smoothed_weights = vector<vector<double> >(N_DIM,vector<double>(N_INTERVALS,0));
    summed_weights = vector<double>(N_DIM,0);
    delta_weights = vector<double>(N_DIM,0);
    average_weight = vector<double>(N_DIM,0);
    std_weight = vector<double>(N_DIM,0);
}
void VEGAS_Map::Reset_Weight()
{
    weights = vector<vector<double> >(N_DIM,vector<double>(N_INTERVALS,0));
    counts = vector<vector<double> >(N_DIM,vector<double>(N_INTERVALS,0));
}

vector<int> VEGAS_Map::Get_Interval_ID(vector<double> y)
{
    vector<int> res;
    for (int i = 0; i < N_DIM; i++)
    {
        res.push_back(floor(y[i]*N_INTERVALS));
    }
    return res;
}
vector<double> VEGAS_Map::Get_Interval_Offset(vector<double> y)
{
    vector<int> ID = Get_Interval_ID(y);
    vector<double> res;
    for (int i = 0; i < N_DIM; i++)
    {
        res.push_back(y[i]*N_INTERVALS - ID[i]);
    }
    return res;
}
vector<double> VEGAS_Map::Get_X(vector<double> y)
{
    vector<int> ID = Get_Interval_ID(y);
    vector<double> offset = Get_Interval_Offset(y);
    vector<double> res;
    for (int i = 0; i < N_DIM; i++)
    {
        int id = ID[i];
        res.push_back(x_edges[i][id] + dx_steps[i][id]*offset[i]);
    }
    return res;
}
double VEGAS_Map::Get_Jac(vector<double> y)
{
    vector<int> ID = Get_Interval_ID(y);
    double jac = 1;
    for (int i = 0; i < N_DIM; i++)
    {
        int id = ID[i];
        jac *= N_INTERVALS*dx_steps[i][id];
    }
    return jac;
}
void VEGAS_Map::Accumulate_Weight(vector<double> y, double f)
{
    // f is the value of integrand!
    vector<int> ID = Get_Interval_ID(y);
    for (int i = 0; i < N_DIM; i++)
    {
        int id = ID[i];
        weights[i][id] += pow(f*Get_Jac(y),2);
        counts[i][id] += 1;
        // cout<<"ID: "<<id<<" weight: "<<weights[i][id]<<" counts: "<<counts[i][id]<<endl;
    }
}
void VEGAS_Map::Smooth_Weight()
{
    // cout<<"Smoothing weight"<<endl;
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        for (int i_inter = 0; i_inter < weights[i_dim].size(); i_inter++)
        {
            if (counts[i_dim][i_inter]!=0)
            {
                weights[i_dim][i_inter] /= counts[i_dim][i_inter];
            }
        }
    }
    // cout<<"Count devided!"<<endl;
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        double d_tmp;
        double d_sum = accumulate(weights[i_dim].begin(),weights[i_dim].end(),0.0);
        summed_weights[i_dim] = 0;
        for (int i = 0; i < N_INTERVALS; i++)
        {
            if (i==0)
            {
                d_tmp = (7.0*weights[i_dim][0]+weights[i_dim][1])/(8.0*d_sum);
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
                d_tmp = (weights[i_dim][N_INTERVALS-2]+7.0*weights[i_dim][N_INTERVALS-1])/(8.0*d_sum);
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
                d_tmp = (weights[i_dim][i-1] + 6.0*weights[i_dim][i] + weights[i_dim][i+1])/(8.0*d_sum);
                if (d_tmp == 0)
                {
                    d_tmp = 0;
                }
                else
                {
                    d_tmp = pow((d_tmp-1.0)/log(d_tmp),alpha);
                }
            }
            smoothed_weights[i_dim][i]=d_tmp;
            summed_weights[i_dim] += d_tmp;
        }
        delta_weights[i_dim] = summed_weights[i_dim]/N_INTERVALS;
    }
}
void VEGAS_Map::Update_Map()
{
    Smooth_Weight();
    // cout<<"Updating the map"<<endl;
    x_edges_last = x_edges;
    dx_steps_last = dx_steps;
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        int current_old_interval = 0;
        int current_new_interval = 1;
        double d_accu = 0;
        while (true)
        {
            d_accu += delta_weights[i_dim];
            while (d_accu > smoothed_weights[i_dim][current_old_interval])
            {
                d_accu -= smoothed_weights[i_dim][current_old_interval];
                current_old_interval++;
            }
            x_edges[i_dim][current_new_interval] = x_edges_last[i_dim][current_old_interval] + d_accu/smoothed_weights[i_dim][current_old_interval]*dx_steps_last[i_dim][current_old_interval];
            dx_steps[i_dim][current_new_interval-1] = x_edges[i_dim][current_new_interval] - x_edges[i_dim][current_new_interval-1];
            current_new_interval++;
            if (current_new_interval >= N_INTERVALS)
            {
                break;
            }
        }
        dx_steps[i_dim][N_INTERVALS-1] = x_edges[i_dim][N_EDGES-1] - x_edges[i_dim][N_EDGES-2];
    }
    Reset_Weight();
}
void VEGAS_Map::Checking_Weight()
{
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        average_weight[i_dim] = 0;
        for (int i = 0; i < weights[i_dim].size(); i++)
        {
            average_weight[i_dim] += weights[i_dim][i];
        }
        average_weight[i_dim] /= weights[i_dim].size();
        for (int i = 0; i < weights[i_dim].size(); i++)
        {
            std_weight[i_dim] += pow(weights[i_dim][i]-average_weight[i_dim],2);
        }
        std_weight[i_dim] = sqrt(std_weight[i_dim]);// /average_weight;
    }
}
void VEGAS_Map::Print_Edges()
{
    cout<<"Grid Map:"<<endl;
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        cout<<"\tx_"<<i_dim<<":";
        for (int i = 0; i < N_EDGES; i++)
        {
            cout<<"\t"<<x_edges[i_dim][i];
        }
        cout<<endl;
        cout<<"\tdx_"<<i_dim<<":";
        for (int i = 0; i < N_INTERVALS; i++)
        {
            cout<<"\t"<<dx_steps[i_dim][i];
        }
        cout<<endl;
    }
}
void VEGAS_Map::Print_Weights()
{
    cout<<"Weights:"<<endl;
    for (int i_dim = 0; i_dim < N_DIM; i_dim++)
    {
        cout<<"\tweight_"<<i_dim<<":";
        for (int i = 0; i < N_INTERVALS; i++)
        {
            cout<<"\t"<<weights[i_dim][i];
        }
        cout<<endl;
        Checking_Weight();
        cout<<"\t\tAverage: "<<average_weight[i_dim]<<"  STD: "<<std_weight[i_dim]<<" std/ave: "<<std_weight[i_dim]/average_weight[i_dim]/N_INTERVALS<<endl;
    }
}
double VEGAS_Map::Checking_Map()
{
    double dx_ave = 1.0/N_INTERVALS;
    double chi2 = 0;
    for (int idim = 0; idim < N_DIM; idim++)
    {
        for (int i = 0; i < N_EDGES; i++)
        {
            chi2 += pow(x_edges[idim][i] - x_edges_last[idim][i],2)/pow(dx_ave,2);
        }
    }
    return chi2/N_DIM/N_EDGES;
}