#include "whale.h"

Whale::Whale(const Input &input, int maxEvaluations, int maxFails)
    : dimension_(input.getDimension()), distances_(input.getDistances()), flow_(input.getFlow()) 
{

    solution_.resize(dimension_);

    // Fix parameters
    maxEvaluations_ = maxEvaluations;
    maxFails_ = maxFails;
    solution_ = randomSolution();
    cost_ = calculateCost(solution_);

    execute();

}

void Whale::execute(){

    

}

int Whale::calculateCost(vector<int> &solution)
{
    int cost = 0;
    for (int i = 0; i < dimension_; ++i)
    {
        for (int j = 0; j < dimension_; ++j)
        {
            // Objetive function
            cost += flow_[i][j] * distances_[solution[i]][solution[j]];
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
        solution[i] = i;
    }

    // Shuffle vector
    random_shuffle (solution.begin(), solution.end());

    return solution;
}