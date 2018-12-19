#ifndef HERD_H__
#define HERD_H__

#include "input.h"
#include "common.h"
#include "whale.h"
#include <deque>

using namespace std;
using namespace std::chrono;

class Herd
{
    public:
        explicit Herd(const Input &input, int evaluations = 100000, int maxFails = 20);

        void execute();

        Whale getBestWhale();

        vector<int> getBestSolution();

        void updateBestSolution();

        int getBestScore();

        void updateParameter(int t);

        Whale getWhale(int indice);

        int compare_function(const void *a,const void *b);

        void swap(vector<int> &solution, int r, int s);

        int moveCost(vector<int> &oldSolution, int r, int s);

        int calculateCost(vector<int> &solution);

        bool checkMove(int r, int s, vector<int> &currentSolution);

        void tabuReset(int &tabuLength);

        //double getTime();

        typedef struct Neighbour
        {
            int r, s, cost;
            bool operator<(const Neighbour &n) const
            {
                return cost < n.cost;
            }
        } Neighbour;

        // Neighbourhood
        vector<Neighbour> neighbourhood_;

        int cost_;
        vector<int>solution_;

        typedef struct
        {
            int i, j, posI, posJ;
        } TabuMove;

        // Tabu Moves
        deque<TabuMove> tabuList_;


    private:

        // Input data
        int dimension_;
        vector< vector<int> > distances_;
        vector< vector<int> > flow_;

        int maxEvaluations_ ;
        int maxFails_ ;

        float a_;
        float r1_;
        float r2_;
        float A_;
        float C_;
        float p_;
        float l_;
        float b_;
        vector<float> D_;
        vector<int> bestSolution_;
        vector<float> D_prime_;
        int bestScore_;
        time_t start_;
        time_t endwait_;
        //Input input_;
        vector<Whale> whale_list_;
       // Whale bestWhale_;
        int taille_;
        int neighbourhoodlength_;
        int tabusize_;

        double time_;

};

#endif