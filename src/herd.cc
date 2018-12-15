#include "herd.h"

Herd::Herd(const Input &input, int maxEvaluations, int maxFails) 
    : dimension_(input.getDimension()), distances_(input.getDistances()), flow_(input.getFlow()), bestWhale_(input) 
{
    //initialization
    //input_ = input;
    maxEvaluations_ = maxEvaluations;
    maxFails_ = maxFails;
    bestSolution_.resize(dimension_);
    start_ = time(NULL);
    time_t seconds = 10; // end loop after this time has elapsed
    endwait_ = start_ =  + seconds;
    taille_ = 10;
    for(int i=0;i<taille_;++i)
    {
        Whale w(input);
        /*whale_list_.add(w);*/
    }

    execute();
    
}

void Herd::execute()
{
while (start_ < endwait_)
    {
        bestWhale_ = getBestWhale();
        for(int i = 0 ; i<taille_ ; ++i)
        {
            start_ = time(NULL);
            printf("loop time is : %s", ctime(&start_));    

            //Update WA parameters (a, A, C, L, p)
            updateParameter();
            
            //D_ = abs(C_*bestWhale_.getSolution() − whale_list_[i].getSolution())
            
            D_ = a_;

            if(p_<0.5) {
                if(abs(A_)<1) {
                    //update the current search agent
                    //agent_t1 = bestWhale_.getSolution() - A_*D_;
                }
                else if (abs(A_)>=1){
                    //select random search agent
                    //update the current search agent
                    //agent_t1 = agent_random - A_*D_;
                }
            }
        
            else if (p_>=0.5) {
                //update the current search agent
                //agent_t1 = D_prime * exp(b_*l) * cos(2*M_PI*l) + bestWhale_.getSolution();
            }

            /*if (recherche_any_repeated_values) {

            }*/
            //largest real value mapping
            //qsort(<arrayname>,<size>,sizeof(<elementsize>),compare_function);

        }
    }

    printf("end time is %s", ctime(&endwait_));
}

Whale Herd::getBestWhale()
{
    int bestCost = whale_list_[0].calculateCost();
    int indice_best_whale = 0;
    for(int i = 1 ; i < taille_ ; ++i)
    {
        if (bestCost>whale_list_[i].calculateCost())
        {
            bestCost = whale_list_[i].calculateCost();
            indice_best_whale = i;
        }
    }
    return whale_list_[indice_best_whale];
}

int Herd::getBestScore()
{
    int bestCost = whale_list_[0].calculateCost();
    for(int i = 1 ; i < taille_ ; ++i)
    {
        if (bestCost>whale_list_[i].calculateCost())
        {
            bestCost = whale_list_[i].calculateCost();
        }
    }
    return bestCost;
}

vector<int> Herd::getBestSolution()
{
    return getBestWhale().getSolution();
}

void Herd::updateBestSolution()
{
    if(getBestScore() < bestScore_)
    {
        bestSolution_ = getBestSolution();
    }
}

void Herd::updateParameter()
{
    a_ =-(float)rand()*2/(float)RAND_MAX; //a = random entre -2 et 0
    r_ = (float)rand()/(float)RAND_MAX;
    A_= 2*a_*r_ - a_;
    C_ = 2*r_;
    p_ = (float)rand()/(float)RAND_MAX;
    l_ =((float(rand()) / float(RAND_MAX)) * (2)) -1;
}

Whale Herd::getWhale(int indice)
{
    return whale_list_[indice];
}

int Herd::compare_function(const void *a,const void *b) {
    int *x = (int *) a;
    int *y = (int *) b;
    return *y - *x;
}