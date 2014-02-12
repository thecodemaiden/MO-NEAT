//
//  MONEAT.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__MONEAT__
#define __SystemGenerator__MONEAT__

#include <vector>
#include <fstream>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>

#include "MNEdge.h"
#include "MNIndividual.h"

#include "NEATCommon.h"
#include "BaseNEAT.h"

class NEATExtendedSpecies : public NEATSpecies {
public:
    std::map<NEATExtendedSpecies *, double> speciesDist;
};

class MONEAT : public BaseNEAT {
protected:
    std::vector<SystemInfo *> population;
    std::vector<NEATExtendedSpecies *>speciesList;

    int stagnantGenerations;
    int maxStagnation;
 
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one
    
    // finds the Pareto fronts of individuals
    void rankSystems();

    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
        
    NEATExtendedSpecies *chooseCompatibleSpecies(NEATExtendedSpecies *species, double maxDist); // fitness proportionate selection
    
public:
    MONEAT(int populationSize, int maxGenerations, int maxStagnation);
    ~MONEAT();
    
    bool tick();

private:
    int nextSpeciesNumber;
    void speciate(); // divide everything into species
    void updateSharedFitnesses();
};


#endif /* defined(__SystemGenerator__MONEAT__) */
