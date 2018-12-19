#include "common.h"
#include "greedy.h"
#include "input.h"
#include "localSearch.h"
#include "tabuSearch.h"
#include "memetic.h"
#include "whale.h"
#include "herd.h"
#include <time.h>

using namespace std;

int main (int argc, char *argv[])
{
   
    srand(time(0));

    vector<int> solution;

    if (argc < 2 || argc > 3)
    {
        cerr << "ERROR: Execute with: qap [data]" << endl;
        return -1;
    }


    // Set seed
    //int seed = 123456;

    // Read file
    Input input(argv[1]);
    if (!input.read())
    {
        cerr << "ERROR: Not file found";
        return -1;
    }
    /*    
    cout << endl;

    // Tabu Search with 100000 evaluations
    //srand(seed);
    TabuSearch ts(input);

    cout << "Tabu Search: " << endl;
    cout << "\tCost: " << ts.getCost() << endl;
    cout << "\tTime: " << ts.getTime() << endl;
    cout << "\tSolution: ";
    vector<int> solution1 =  ts.getSolution();
    for (int i = 0; i < solution1.size(); ++i)
    {
        cout << solution1[i] << " ";
    } cout << endl;

    cout << endl;*/
/*
    // Stationary Genetic Algorythm with position based crossings with 100000 evaluations
    //srand(seed);
    Memetic m(input, 1000000, 20, 0.2);

    cout << "Memetic: " << endl;
    cout << "\tCost: " << m.getCost() << endl;
    cout << "\tTime: " << m.getTime() << endl;
    cout << "\tSolution: ";
    solution1 =  m.getSolution();
    for (int i = 0; i < solution.size(); ++i)
    {
        cout << solution[i] << " ";
    } cout << endl;

    cout << endl;*/
 /*
    Whale w(input);

    cout << "Whale: " << endl;
    cout << "\tCost: " << w.getCost() << endl;
    //cout << "\tTime: " << w.getTime() << endl;
    cout << "\tSolution: ";
    solution =  w.getSolution();
    for (int i = 0; i < solution.size(); ++i)
    {
        cout << solution[i] << " ";
    } cout << endl;

*/
double temps_total = 0;
for(int i= 0 ; i<10 ; ++i)
{
    cout<<i+1<<endl;
    Herd h(input);
    /*cout<<"Time : ";
    cout<<h.getTime()<<endl;
    cout<<endl;
    temps_total+=h.getTime();*/
}
//cout<<"temps total = "<<temps_total<<endl;

    return 0;
}