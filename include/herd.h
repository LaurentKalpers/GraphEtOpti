#ifndef HERD_H__
#define HERD_H__

#include "input.h"
#include "common.h"
#include "whale.h"

using namespace std;

class herd
{
    public:
        explicit Herd(Input input, int evaluations = 100000, int maxFails = 20);

        void execute();

        Whale getBestWhale();

        vector<int> getBestSolution();

        void updateBestSolution();

        int getBestScore();

        void updateParameter();

        Whale getWhale(int indice);

        int compare_function(const void *a,const void *b);


    private:
        float a_;
        float r_;
        float A_;
        float C_;
        float p_;
        float l_;
        float D_;
        vector<int> solution;
        int bestScore;
        time_t start_;
        time_t endwait_;
        Input input_;
        vector<Whale> whale_list_;
        Whale bestWhale_;
        int taille_;

};

#endif