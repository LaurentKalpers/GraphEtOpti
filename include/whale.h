#ifndef WHALE_H__
#define WHALE_H__

#include "input.h"
#include "common.h"

using namespace std;

class Whale
{
    public: 

        explicit Whale(const Input &input, int evaluations = 100000, int maxFails = 20);

        // Calculate cost using objetive function
        int calculateCost(vector<int> &solution);

        void execute();
        

        // Random solution
        vector<int> randomSolution();

        vector<int> getSolution();
        int getCost();
        //double getTime();

    private:
        // Input data
        int dimension_;
        vector< vector<int> > distances_;
        vector< vector<int> > flow_;

        // Solution data
        vector<int> solution_;
        int cost_;
        double time_;

        // Parameters
        int maxEvaluations_;
        int maxFails_;
};

#endif