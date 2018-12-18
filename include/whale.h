#ifndef WHALE_H__
#define WHALE_H__

#include "input.h"
#include "common.h"

using namespace std;

class Whale
{
    public: 

        //  explicit Whale::Whale(const Whale &whale);
        explicit Whale( int dimension, vector< vector<int> > distances, vector< vector<int> > flow);

        // Calculate cost using objetive function
        int calculateCost();//vector<int> &solution);

        void execute();
        

        // Random solution
        vector<int> randomSolution();

        vector<int> getSolution();
        int getCost();
        double getTime();

        void setSolution(vector<int> v);

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