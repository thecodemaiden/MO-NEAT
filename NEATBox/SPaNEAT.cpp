//
//  SPaNEAT.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/9/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "SPaNEAT.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>
#include <assert.h>
#include <map>
#include <iostream>

#include "NBUtils.h"

#define USE_NODE_DIFF 1
#define SEPARATE_ARCHIVE 1

SPaNEAT::SPaNEAT(long populationSize, long archiveSize, long maxGenerations)
:BaseNEAT(populationSize, maxGenerations)
{}

SPaNEAT::~SPaNEAT()
{
    //empty the archive
#if SEPARATE_ARCHIVE
    for (std::vector<SPSystemInfo *>::iterator it = archive.begin(); it != archive.end(); it++) {
        delete *it;
    }
    
#endif
    //empty the population
    for (std::vector<SPSystemInfo *>::iterator it = population.begin(); it != population.end(); it++) {
        delete *it;
    }
 
    delete origin;
}

// return true if sys1 comes before sys2
// smaller rank fitness is better
static bool compareIndividuals(SystemInfo *sys1, SystemInfo *sys2)
{
    return sys1->rankFitness < sys2->rankFitness;
}

static bool compareInnovationNumbers(const InnovationInfo *a1, const InnovationInfo *a2)
{
    return a1->innovationNumber < a2->innovationNumber;
}



void SPaNEAT::spawnNextGeneration()
{
    std::vector<SPSystemInfo *>newGeneration;
    
    // create population children through recombination and mutation of archive members
    while (newGeneration.size() < populationSize-archiveSize) {
        const int tournamentSize = 3;
        // select parent one
        SPSystemInfo *p1 = archive[uniformlyDistributed(archiveSize)];
        for (int i=0; i<tournamentSize-1; i++) {
             SPSystemInfo *contender =  archive[uniformlyDistributed(archiveSize)];
            if (contender->rankFitness < p1->rankFitness)
                p1 = contender;
        }
        
        SPSystemInfo *p2 = archive[uniformlyDistributed(archiveSize)];
        for (int i=0; i<tournamentSize-1; i++) {
            SPSystemInfo *contender =  archive[uniformlyDistributed(archiveSize)];
            if (contender->rankFitness < p2->rankFitness)
                p2 = contender;
        }
        
        MNIndividual *i = combineSystems(p1, p2);
        mutateSystem(i);

        SPSystemInfo *child = new SPSystemInfo(i);
        newGeneration.push_back(child);
    }
    
    // add the archive too (for now)
#if SEPARATE_ARCHIVE
    population = newGeneration;
#else
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());
#endif
}


void SPaNEAT::prepareInitialPopulation()
{
    // we need to index the genome of the 'base' system
    origin = createInitialIndividual();
    std::vector<InnovationInfo *> baseGenome = orderedConnectionGenome(origin->connectionGenome());
    
    // clean out the allocated InnovationInfos
    for (std::vector<InnovationInfo *>::iterator dit = baseGenome.begin(); dit != baseGenome.end();) {
        delete *dit;
        dit = baseGenome.erase(dit);
    }
    
    
    // mutate the initial system to get an initial population
    while (population.size() < populationSize) {
        MNIndividual *newSystem = createInitialIndividual();
        mutateSystem(newSystem);
        SPSystemInfo *info = new SPSystemInfo(newSystem);
        population.push_back(info);
    }
}

// does s1 dominate s2?
// i.e. all s1 fitness values are greater than s2?
static int domination(SPSystemInfo *s1, SPSystemInfo *s2)
{
    bool firstDominates = true;
    bool secondDominates = true;
    
    std::vector<double>::iterator fitnessIter1 = s1->fitnesses.begin();
    std::vector<double>::iterator fitnessIter2 = s2->fitnesses.begin();
    
    while (fitnessIter1 != s1->fitnesses.end()) {
        if (*fitnessIter1 > *fitnessIter2)
            firstDominates = false;
        
        if (*fitnessIter2 > *fitnessIter1)
            secondDominates = false;
        fitnessIter2++;
        fitnessIter1++;
    }
    
    if (firstDominates == secondDominates)
        return 0;
    
    if (firstDominates)
        return 1;
    
    if (secondDominates)
        return -1;
    
    return 0;
}

// using SPEA2
void SPaNEAT::rankSystems()
{
    std::vector<SPSystemInfo *>::iterator popIter;
    
    std::map<SPSystemInfo *, std::vector<SPSystemInfo *> > dominationMap;
    
    int k = sqrt(populationSize+archiveSize);
#if SEPARATE_ARCHIVE
    std::vector<SPSystemInfo *> toRank = std::vector<SPSystemInfo *>(population);
    toRank.insert(toRank.end(), archive.begin(), archive.end());
#else
    std::vector<SPSystemInfo *> toRank = population;
#endif
    // we're assigning a strength number based on # of dominated individuals
    // domination map and count take  opposite meanings from NSGA-II
    
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        std::vector<SPSystemInfo *> dominators;
        SPSystemInfo *sys = *popIter;
        
        sys->dominationCount = 0;
        sys->distances.clear();
        
        typename std::vector<SPSystemInfo *>::iterator innerIter;
        for (innerIter = toRank.begin(); innerIter != toRank.end(); innerIter++) {
            if (*innerIter == *popIter)
                continue;
            SPSystemInfo *otherSys = *innerIter;
            
            int dominationScore = domination( otherSys, sys);
            
            if (dominationScore == 1) {
                // *innerIter dominates
                dominators.push_back(otherSys);
            }
            if (dominationScore == -1) {
                // *popIter dominates
                sys->dominationCount += 1;
            }
            sys->distances.push_back(genomeDistance(sys->individual, otherSys->individual));
        }
        dominationMap[sys] = dominators;
        
        // the density calculation can come here, since we have all the distances any way
        std::sort(sys->distances.begin(), sys->distances.end());
        double kDist = sys->distances[k];
        sys->rankFitness = 1.0/(kDist + 2);
    }
    
    
    // now assign raw fitness based on the strengths of dominators, plus population neighborhood density
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        SPSystemInfo *sys = *popIter;
        std::vector<SPSystemInfo *> dominators = dominationMap[sys];
        std::vector<SPSystemInfo *>::iterator it;
        for (it = dominators.begin(); it != dominators.end(); it++) {
            sys->rankFitness += (*it)->dominationCount;
        }
    }
}

bool SPaNEAT::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
        first_run = true;
    }
    
    std::vector<SPSystemInfo *>::iterator popIter;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        SPSystemInfo *sys = *popIter;
        sys->fitnesses.clear();
        
        std::vector<EvaluationFunction>::iterator funcIter;
        for (funcIter =  evaluationFunctions.begin(); funcIter != evaluationFunctions.end(); funcIter++) {
            double val = (*funcIter)(sys->individual);
            sys->fitnesses.push_back(val);
        }
    }
    
    rankSystems();
    
    bool last_run = (generations >= maxGenerations) ;
    
    std::vector<SPSystemInfo *> newArchive;
    
#if SEPARATE_ARCHIVE
    std::vector<SPSystemInfo *> toRank = std::vector<SPSystemInfo *>(population);
    toRank.insert(toRank.end(), archive.begin(), archive.end());
#else
    std::vector<SPSystemInfo *> toRank = population;
#endif
    
    std::sort(toRank.begin(), toRank.end(), compareIndividuals);
    
    std::vector<SPSystemInfo *>::iterator lastArchived = population.begin();
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        if ((*popIter)->rankFitness < 1.0) {
            newArchive.push_back(new SPSystemInfo(**popIter));
            lastArchived = popIter;
        }
    }
    
    std::cout << newArchive.size() << " non dominated solutions (" << generations << ").\n";

    
    long archiveExcess = newArchive.size() - archiveSize;
    
    if (archiveExcess < 0) {
        // add in the next best -archiveExcess solutions
        for (popIter = lastArchived+1; popIter != lastArchived - archiveExcess +1; popIter++) {
            newArchive.push_back(new SPSystemInfo(**popIter));
        }
    }
    
    if (archiveExcess > 0) {
        // remove by crowding distance...
        std::sort(newArchive.begin(), newArchive.end(), [](const SPSystemInfo * elem1, const SPSystemInfo *elem2){return elem1->distances < elem2->distances;});
        
        popIter = newArchive.begin();
        for (long i=0; i<archiveExcess; i++) {
            delete *popIter;
            popIter = newArchive.erase(popIter);
        }
    }

    generations++;

    // delete everything in the old population
    // XXX: this includes the old archive (for now)
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        delete  *popIter;
    }
    population.clear();
    
#if SEPARATE_ARCHIVE
    for (popIter = archive.begin(); popIter != archive.end(); popIter++) {
        delete  *popIter;
    }
#endif
    
    archive = newArchive;
    
    if (!last_run) {
        spawnNextGeneration();
    }
    
    return  last_run;
}

std::vector<SystemInfo *> SPaNEAT::optimalSolutions()
{
    // some suboptimal solutions might be in here: make sure we only get the nondominated ones
        std::vector<SystemInfo *>converted;
    
    std::vector<SPSystemInfo *>::iterator it = std::find_if(archive.begin(), archive.end(), [](const SystemInfo * elem){return elem->rankFitness < 1.0;});
    for (; it != archive.end(); it = std::find_if(++it, archive.end(), [](const SystemInfo * elem){return elem->rankFitness < 1.0;})) {
        converted.push_back(*it);
    }

    return converted;
}