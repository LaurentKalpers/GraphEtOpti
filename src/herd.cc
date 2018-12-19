#include "herd.h"
#include <math.h> 
#include <time.h>
#include <iostream>
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

            if(calculateCost(solutiontmp)< whale_list_[i].getCost()){
            whale_list_[i].setSolution(solutiontmp);
            }

            cout << "\t solution whale :  ";
            for (int j = 0; j < dimension_; ++j)
            {
                cout << whale_list_[i].getSolution()[j] << " ";
            } cout << endl;
            cout<<whale_list_[i].getCost()<<endl;

            agent_t1.clear();
            D_.clear();
            D_prime_.clear();
           // D_.resize(0);
            //D_prime_.resize(0);
            
            

            //TABU SEARCH
            // Empty the neighborhoud
				int neighbourhoodSize_;
				neighbourhoodSize_ = 30;
                vector<int>currentSolution=whale_list_[i].getSolution();
                cost_=whale_list_[i].getCost();
                int currentCost=cost_;
                bool improvement=false;
                int tabuLength = dimension_ / 5;

			// Generate the Neighbourhood
			for (int i = 0; i < neighbourhoodSize_; ++i)
			{
				// Generate random r and s
				int r = rand() % dimension_;
				int s = rand() % dimension_;



				// Check the difference
				int diff = moveCost(currentSolution, r, s);

				// Add to the neighbourhood
				Neighbour neighbor;
				neighbor.r = r;
				neighbor.s = s;
				//neighbor.cost = diff + currentCost;
				neighbor.cost = diff + whale_list_[i].getCost();

				neighbourhood_.push_back(neighbor);
			}

			// Sort the neighbourhood
			sort(neighbourhood_.begin(), neighbourhood_.end());

			// See the neighbourhood
			bool selected = false;
			int z = 0;

			// Check if the move is not Tabu
			selected = checkMove(neighbourhood_[z].r, neighbourhood_[z].s, currentSolution);

			// Check if is better than our solution
			int nCost = whale_list_[i].getCost() + moveCost(currentSolution, neighbourhood_[z].r, neighbourhood_[z].s);
			if (!selected && nCost < cost_)
			{
				selected = true;
			}
			z++;

			while (!selected && z < neighbourhoodSize_)
			{
				// Check if the move is not Tabu
				selected = checkMove(neighbourhood_[z].r, neighbourhood_[z].s, currentSolution);
				z++;
			}	

			// If none of the neighbourhood is selected, generate new neighbourhood
			if (!selected)
			{
				continue;
			}

			int r = neighbourhood_[z - 1].r;
			int s = neighbourhood_[z - 1].s;

            cout << "\t Avant swap :  ";
            for (int j = 0; j < dimension_; ++j)
            {
                cout << currentSolution[j] << " ";
            } cout << endl;
            cout<<calculateCost(currentSolution)<<endl;

			// Jump to something
			swap(currentSolution, r, s);
            cout << "\t Après swap :  ";
            for (int j = 0; j < dimension_; ++j)
            {
                cout << currentSolution[j] << " ";
            } cout << endl;
            cout<<calculateCost(currentSolution)<<endl;

            whale_list_[i].setSolution(currentSolution);
            //cout<<"Après set solution : "<<whale_list_[i].getSolution()<<"\n";
			currentCost = neighbourhood_[z - 1].cost;
            //cout<<"Après set solution : "<<whale_list_[i].getCost()<<"\n";

			// Update frecuency matrix
			//frec_[r][currentSolution[r]]++;
			//frec_[s][currentSolution[s]]++;

			if (currentCost < cost_)
			{
				improvement = true;
				solution_ = currentSolution;
				cost_ = whale_list_[i].getCost();
			}

			// Create tabu move
			TabuMove tabuMove;
			tabuMove.i = currentSolution[neighbourhood_[z - 1].r];
			tabuMove.j = currentSolution[neighbourhood_[z - 1].s];
			tabuMove.posI = neighbourhood_[z - 1].r;
			tabuMove.posJ = neighbourhood_[z - 1].s;

			// Add tabu move
			tabuList_.push_back(tabuMove);

			//Size of the deque
			if (tabuList_.size() > tabuLength)
			{
				tabuList_.pop_front();
			}

            neighbourhood_.clear();






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
    cout << "\tmeilleures solution tabu :  ";
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
void Herd::swap(vector<int> &solution, int r, int s)
{
	int t;
	// Simple swap
	t = solution[r];
	solution[r] = solution[s];
	solution[s] = t;
}

int Herd::moveCost(vector<int> &oldSolution, int r, int s)
{
	int cost = 0;
	for (int k = 0; k < dimension_; ++k)
	{
		if (k != r && k != s)
		{
			cost +=
				flow_[r][k] * (distances_[oldSolution[s]][oldSolution[k]] - distances_[oldSolution[r]][oldSolution[k]]) +
				flow_[s][k] * (distances_[oldSolution[r]][oldSolution[k]] - distances_[oldSolution[s]][oldSolution[k]]) +
				flow_[k][r] * (distances_[oldSolution[k]][oldSolution[s]] - distances_[oldSolution[k]][oldSolution[r]]) +
				flow_[k][s] * (distances_[oldSolution[k]][oldSolution[r]] - distances_[oldSolution[k]][oldSolution[s]]);
		}
	}

	return cost;
}

int Herd::calculateCost(vector<int> &solution)
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

// Search in the Tabu Moves if the move r,s is allowed with the vector current solution
bool Herd::checkMove(int r, int s, vector<int> &currentSolution)
{
	for (int n = 0; n < tabuList_.size(); ++n)
	{
		if ((tabuList_[n].i == r && tabuList_[n].posI == currentSolution[r]))
		{
			return false;
		}

		if ((tabuList_[n].i == s && tabuList_[n].posI == currentSolution[s]))
		{
			return false;
		}

		if ((tabuList_[n].j == r && tabuList_[n].posJ == currentSolution[s]))
		{
			return false;
		}

		if ((tabuList_[n].j == s && tabuList_[n].posJ == currentSolution[s]))
		{
			return false;
		}
	}

	return true;
}

void Herd::tabuReset(int &tabuLength)
{
	// Clear Tabu List and change tabu size
	tabuList_.clear();

	float random = ((double)rand() / (RAND_MAX));
	if (random < 0.5)
	{
		tabuLength = tabuLength + (0.5 * tabuLength);
	}
	else
	{
		tabuLength *= random - (0.5 * tabuLength);
	}
}