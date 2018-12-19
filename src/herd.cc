#include "herd.h"
#include <math.h> 
#include <time.h>
# define M_PI           3.14159265358979323846

Herd::Herd(const Input &input, int maxEvaluations, int maxFails) 
    : dimension_(input.getDimension()), distances_(input.getDistances()), flow_(input.getFlow())
{
    //initialization
    //input_ = input;
    maxEvaluations_ = maxEvaluations;
    maxFails_ = maxFails;
    bestSolution_.resize(dimension_);
    time_t start_ = time(NULL);
    time_t seconds = 10000000; // end loop after this time has elapsed
    endwait_ = start_ =  + seconds;

    b_ = 1;
    taille_ = 20;

    for(int i=0;i<taille_;++i)
    {
        Whale w(dimension_, distances_, flow_);
        cout<<w.getCost()<<endl;
        for (int j = 0; j < dimension_; ++j)
        {
        cout << w.getSolution()[j] << " ";
        } cout << endl;
        whale_list_.push_back(w);
    }

    bestScore_ = getBestScore();
    cout << "\tCost: " << bestScore_ << endl;      
    cout << "\tSolution: ";
    for (int i = 0; i < dimension_; ++i)
    {
        cout << getBestWhale().getSolution()[i] << " ";
    } cout << endl;

    execute();
    bestScore_ = getBestScore();
    cout << "\tCost: " << bestScore_ << endl;      
    cout << "\tSolution: ";
    for (int i = 0; i < dimension_; ++i)
    {
        cout << getBestWhale().getSolution()[i] << " ";
    } cout << endl;
}

void Herd::execute()
{
    int randNum;
    int t=0;
    Whale bestWhale_ = getBestWhale();
    vector<int> solutiontmp;
    //solutiontmp.reserve(dimension_);
    vector<float> agent_t1;
    //while (start_ < endwait_)
    Whale bestbestwhale(dimension_, distances_, flow_);
    while(t<5001)
    {
        cout<<t<<endl;
        for(int i = 0 ; i<taille_ ; ++i)
        {
            solutiontmp.resize(dimension_);

            start_ = time(NULL);

            //Update WA parameters (a, A, C, L, p)
            updateParameter(t);   

            for (int j=0 ; j<dimension_ ; ++j)
            {
                D_.push_back(abs(C_*bestWhale_.getSolution()[j] - whale_list_[i].getSolution()[j]));
                D_prime_.push_back(abs(bestWhale_.getSolution()[j] - whale_list_[i].getSolution()[j]));
            }

            //cout<<"update"<<endl;
            
            if(p_<0.5) 
            {    
                if(abs(A_)<1) 
                {
                    //update the current search agent
                    //agent_t1 = bestWhale_.getSolution() - A_*D_;
                    for (int j=0 ; j<dimension_ ; ++j)
                    {
                        agent_t1.push_back(bestWhale_.getSolution()[j] - A_*D_[j]);
                    }
                }
                else if (abs(A_)>=1)
                {
                    //select random search agent
                    //update the current search agent
                    //agent_t1 = agent_random - A_*D_;

                    //Whale w(dimension_, distances_, flow_);
                    randNum = rand()%(taille_);// + 1) ;
                    for (int j=0 ; j<dimension_ ; ++j)
                    {
                        //agent_t1.push_back(w.getSolution()[j] -  A_*D_[j]);
                        agent_t1.push_back(whale_list_[randNum].getSolution()[j] -  A_*D_[j]);
                    }
                }
            }
            
            else if (p_>=0.5) 
            {
                //update the current search agent
                //agent_t1 = D_prime * exp(b_*l) * cos(2*M_PI*l) + bestWhale_.getSolution();
                for (int j=0 ; j<dimension_ ; ++j)
                {
                    agent_t1.push_back(D_prime_[j] * exp(b_*l_) * cos(2*M_PI*l_) + bestWhale_.getSolution()[j]);
                }
            }

            /*if (recherche_any_repeated_values) {}*/

            //largest real value mapping
            //qsort(<arrayname>,<size>,sizeof(<elementsize>),compare_function);
            //qsort(agent_t1, sizeof(agent_t1),sizeof(agent_t1),&Herd::compare_function);

            for(int j = 0; j<dimension_;++j)
            {         
                int maxElementIndex = max_element(agent_t1.begin(),agent_t1.end()) - agent_t1.begin();
                agent_t1[maxElementIndex] = -100000000000000 ;
                solutiontmp[maxElementIndex] = dimension_-j-1;   
                //cout<<"real mapping"<<endl;       
            }


            whale_list_[i].setSolution(solutiontmp);

            agent_t1.clear();
            D_.clear();
            D_prime_.clear();
           // D_.resize(0);
            //D_prime_.resize(0);
            bestWhale_ = getBestWhale();
            if(bestWhale_.calculateCost()<bestbestwhale.calculateCost()){
                bestbestwhale = bestWhale_;
            }
            else {
                bestWhale_ = bestbestwhale;
            }

        }
        t=t+1;
    }
            cout << "\tmeilleures solution :  ";
            for (int j = 0; j < dimension_; ++j)
            {
                cout << bestbestwhale.getSolution()[j] << " ";
            } cout << endl;
            cout<<bestbestwhale.getCost()<<endl;
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

void Herd::updateParameter(int t)
{ 
    a_ =2-t/2500; 
    r1_ = (float)rand()/(float)RAND_MAX;
    r2_ = (float)rand()/(float)RAND_MAX;
    A_= 2*a_*r1_ - a_;
    C_ = 2*r2_;
    p_ = (float)rand()/(float)RAND_MAX;
    l_ =((float(rand()) / float(RAND_MAX)) * (2)) -1;
}

Whale Herd::getWhale(int indice)
{
    return whale_list_[indice];
}

int Herd::compare_function(const void *a,const void *b) 
{
    int *x = (int *) a;
    int *y = (int *) b;
    return *y - *x;
}