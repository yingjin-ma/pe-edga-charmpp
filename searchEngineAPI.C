#ifndef __SEARCHENGINEAPI__
#define __SEARCHENGINEAPI__

#include "searchEngine.h"
//#include "parameters.h"

/*   framework for search engine */

extern int initial_grainsize;
extern int numQueens;
//extern CalcParams Params; 

class GAStateBase : public StateBase
{
public:
    
};

void createInitialChildren(Solver *solver)
{

}


inline void createChildren( StateBase *_base , Solver* solver, bool parallel)
{

}

    inline double cost( )
    {
        return 0;
    }

    double heuristic( )
    {
        return 0;
    }

    double bound( int &l )
    {
        return 0;
    }

    //Search Engine Option
    int parallelLevel()
    {
        return initial_grainsize;
    }

    int searchDepthLimit()
    {
        return 1;
    }

    int minimumLevel()
    {
        return 1;
    }

    int maximumLevel()
    {
        return  numQueens;
    }
    inline void searchDepthChangeNotify( int ) {}


    SE_Register(GAStateBase, createInitialChildren, createChildren, parallelLevel, searchDepthLimit);

#endif
