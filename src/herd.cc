#include "herd.h"

Herd::Herd(Input input, int maxEvaluations, int maxFails)
{
    //initialization
    input_ = input;
    time_t start_ = time(NULL);
    time_t seconds = 10; // end loop after this time has elapsed
    endwait_ = start + seconds;
    taille_ = 10;
    for(int i=0;i<taille_;++i)
    {
        Whale w(input_);
        whale_list_.add(w);
    }
    
}

void Herd::execute()
{
while (start < endwait)
    {
        bestWhale_ = getBestWhale;
        for(int i = 0 ; i<taille_ ; ++i)
        {
            /* Do stuff while waiting */
            start = time(NULL);
            printf("loop time is : %s", ctime(&start));    

            //Update WA parameters (a, A, C, L, p)
            updateParameter();
            
            D_ = abs(C_*bestWhale_.getSolution() − whale_list_[i].getSolution())

            if(p<0.5) {
                if(abs(A)<1) {
                    //update the current search agent
                    agent_t1 = bestWhale_.getSolution() - A_*D_;
                }
                else if (abs(A)>=1){
                    //select random search agent
                    //update the current search agent
                    agent_t1 = agent_random - A_*D_;
                }
            }
            else if (p>=0.5) {
                //update the current search agent
                agent_t1 = D_prime * exp(b_*l) * cos(2*M_PI*l) + bestWhale_.getSolution();
            }

            if (recherche_any_repeated_values) {

            }
            //largest real value mapping
            qsort(<arrayname>,<size>,sizeof(<elementsize>),compare_function);

        }
    }

    printf("end time is %s", ctime(&endwait));
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
        bestScore_ = getBestSolution();
    }
}

void Herd::updateParameter()
{
    a_ =-(float)rand()*2/(float)RAND_MAX; //a = random entre -2 et 0
    r_ = (float)rand()/(float)RAND_MAX;
    A_= 2*a*r - a;
    C_ = 2*r;
    p_ = (float)rand()/(float)RAND_MAX;
    l_ =((float(rand()) / float(RAND_MAX)) * (2)) -1;
}

Whale Herd::getWhale(int indice)
{
    return whale_list_[indice];
}

int compare_function(const void *a,const void *b) {
    int *x = (int *) a;
    int *y = (int *) b;
    return *y - *x;
}