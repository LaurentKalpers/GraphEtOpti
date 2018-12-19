#include "whale.h"

Whale::Whale(int dimension, vector< vector<int> > distances,vector< vector<int> > flow)
{

    dimension_ = dimension;
    distances_ = distances;
    flow_ = flow;

    solution_.resize(dimension_);

    // Fix parameters
    solution_ = randomSolution();
    cost_ = calculateCost();

    execute();

}


void Whale::setSolution(vector<int> solution)
{
    solution.resize(dimension_);
    solution_.resize(dimension_);
    solution_ = solution;
}

void Whale::execute(){
}

int Whale::calculateCost()//vector<int> &solution)
{
    int cost = 0;
    for (int i = 0; i < dimension_; ++i)
    {
        for (int j = 0; j < dimension_; ++j)
        {
            // Objetive function
            cost += flow_[i][j] * distances_[solution_[i]][solution_[j]];
        }
    }

    return cost;
}

vector<int> Whale::getSolution()
{
    return solution_;
}

int Whale::getCost()
{
    cost_ = calculateCost();
    return cost_;
}

double Whale::getTime()
{
    return time_;
}

vector<int> Whale::randomSolution()
{
    vector<int> solution;
    solution.resize(dimension_);

    // Initialize
    for (int i = 0; i < dimension_; ++i)
    {
        solution[i] = i; //mettre un +1 dans l affichage
    }

    // Shuffle vector
    random_shuffle (solution.begin(), solution.end());

    return solution;
}
