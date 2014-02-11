
//
//  MONEAT.cpp
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#include "MONEAT.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>
#include <assert.h>
#include <map>
#include <iostream>

#include "NBUtils.h"

#define RESTRICT_SPECIES 0
#define ROULETTE_SELECT 1
#define INTERSPECIES 1

#define TIGHT_CLUSTERS 0

#define USE_NODE_DIFF 1

MONEAT::MONEAT(int populationSize, int maxGenerations, int maxStagnation)
:BaseNEAT(populationSize, maxGenerations), nextSpeciesNumber(1),
stagnantGenerations(0)
{
}

void MONEAT::prepareInitialPopulation()
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
        SystemInfo *info = new SystemInfo(newSystem);
        population.push_back(info);
    }
}

MONEAT::~MONEAT()
{
    //empty the species lists
    std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        delete *speciesIterator;
        speciesIterator = speciesList.erase(speciesIterator);
    }
    
    // empty the population
    std::vector<SystemInfo *>::iterator it = population.begin();
    while (it != population.end()) {
        delete *it;
        it = population.erase(it);
    }
    delete origin;
}

// return true if sys1 comes before sys2
static bool compareIndividuals(SystemInfo *sys1, SystemInfo *sys2)
{
    return sys1->rankFitness > sys2->rankFitness;
}

static bool distanceCompare(const std::pair<NEATExtendedSpecies *, double> &e1, const std::pair<NEATExtendedSpecies *, double> &e2)
{
    return e1.second < e2.second;
}

void MONEAT::spawnNextGeneration()
{
 
#if INTERSPECIES
    double maxSpeciesDist = -INFINITY;
    for (std::vector<NEATExtendedSpecies *>::iterator speciesIt = speciesList.begin(); speciesIt != speciesList.end(); speciesIt++) {
        double maxD = std::max_element((*speciesIt)->speciesDist.begin(), (*speciesIt)->speciesDist.end(), distanceCompare)->second;
        if (maxD > maxSpeciesDist)
            maxSpeciesDist = maxD;
    }
#endif
    
    // create population children through recombination
    std::vector<SystemInfo *>newGeneration;
    
    
    std::vector<NEATExtendedSpecies *>::iterator speciesIt = speciesList.begin();
    while(speciesIt != speciesList.end()) {
        // selectionnnnn
        NEATExtendedSpecies *breedingSpecies = *speciesIt;
        std::vector<SystemInfo *>newMembers;
        
        std::vector<SystemInfo *>individuals = breedingSpecies->members;
        
        while (newMembers.size() < individuals.size()) {
            
            SystemInfo *p1;
            SystemInfo *p2;
            
            MNIndividual *child;
            
#if ROULETTE_SELECT
            //stochastic acceptance
            bool accepted = false;
            while (!accepted) {
                p1 = individuals[uniformlyDistributed(individuals.size())];
                accepted = (uniformProbability() < p1->rankFitness);
            }
            
#if INTERSPECIES
            std::vector<SystemInfo *> individuals2 = chooseCompatibleSpecies(breedingSpecies, maxSpeciesDist)->members;
#else
            std::vector<SystemInfo *>individuals2 = individuals;
#endif
            
            accepted = false;
            while (!accepted) {
                p2 = individuals2[uniformlyDistributed(individuals2.size())];
                accepted = (uniformProbability() < p2->rankFitness);
            }
            
#else
            const int tournamentSize = 3;
  
#if INTERSPECIES
            std::vector<SystemInfo *> individuals2 = chooseCompatibleSpecies(breedingSpecies, maxSpeciesDist)->members;
#else
            std::vector<SystemInfo *>individuals2 = individuals;
#endif
            
            // tournament
            p1 = (individuals[uniformlyDistributed(individuals.size())]);
            p2 = (individuals2[uniformlyDistributed(individuals2.size())]);
            
            
            if (individuals.size() > tournamentSize) {
                for (int i=0; i<tournamentSize; i++) {
                    SystemInfo *contender = individuals[uniformlyDistributed(individuals.size())];
                    if (contender->rankFitness > p1->rankFitness)
                        p1 = contender;
                }
                
                for (int i=0; i<tournamentSize; i++) {
                    SystemInfo *contender = individuals2[uniformlyDistributed(individuals2.size())];
                    if (contender->rankFitness > p2->rankFitness)
                        p2 = contender;
                }
            }
#endif
            child = combineSystems(p1, p2);
            // mutate the new individual (maybe)
            mutateSystem(child);
            
            SystemInfo *i = new SystemInfo(child);
            newGeneration.push_back(i);
            newMembers.push_back(i);
        }
        (*speciesIt)->members.insert((*speciesIt)->members.end(), newMembers.begin(), newMembers.end());
        speciesIt++;
    }
    // add the original population too
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());
    
}

NEATExtendedSpecies *MONEAT::chooseCompatibleSpecies(NEATExtendedSpecies *species, double maxDistance){
    NEATExtendedSpecies *chosen = NULL;
    bool accepted = false;
    while (!accepted) {
        chosen = speciesList[uniformlyDistributed(speciesList.size())];
        double d = species->speciesDist[chosen];
        accepted = (uniformProbability() < ((maxDistance - d)/maxDistance));
    }
    
    return chosen;
}


void MONEAT::speciate()
{
    std::vector<SystemInfo *>::iterator populationIter = population.begin();
    while (populationIter != population.end()) {
        
        // assign to species
        bool added= false;
#if TIGHT_CLUSTERS
        double smallestDiff = INFINITY;
        std::vector<NEATExtendedSpecies *>::iterator chosenSpecies = speciesList.end();
#endif
        
        std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            double  distance = fabs(genomeDistance((*populationIter)->individual, (*speciesIterator)->representative));
            if (distance < d_threshold) {
#if TIGHT_CLUSTERS
                if (distance < smallestDiff) {
                    smallestDiff = distance;
                    chosenSpecies = speciesIterator;
                }
#else
                added = true;
                (*speciesIterator)->members.push_back((*populationIter));
                break;
#endif
            }
            speciesIterator++;
        }
#if TIGHT_CLUSTERS
        if (chosenSpecies != speciesList.end()) {
            added = true;
            (*chosenSpecies)->members.push_back((*populationIter));
        }
#endif
        
        if (!added) {
            // new species!!
            NEATExtendedSpecies *s = new NEATExtendedSpecies();
            s->representative = (*populationIter)->individual->clone();
            s->members.push_back(*populationIter);
            speciesList.push_back(s);
            s->speciesNumber = nextSpeciesNumber++;
        }
        
        populationIter++;
    }
}

void MONEAT::updateSharedFitnesses()
{
    long nSpecies = speciesList.size();
    fprintf(stderr, "%ld species\n", nSpecies);
    
    // for each species, copy a random current member to be the representative for the next generation, and adjust the fitnesses for sharing
    // kill off a species with no members
    std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        // pick a new rep
        std::vector<SystemInfo *> members = (*speciesIterator)->members;
        if (members.size() == 0) {
            NEATExtendedSpecies *extinctSpecies = *speciesIterator;
            delete extinctSpecies;
            speciesIterator = speciesList.erase(speciesIterator);
        } else {
            long index = uniformlyDistributed(members.size());
            MNIndividual *rep = (members[index]->individual);
            delete (*speciesIterator)->representative;
            (*speciesIterator)->representative = rep->clone();
            
            // fprintf(stderr, "\tSpecies %d: %ld members\n", ++i, members.size());
            
            // update the fitnesses
            double  totalFitness = 0.0;
            typename std::vector<SystemInfo *>::iterator memberIterator = members.begin();
            while (memberIterator != members.end()) {
                totalFitness += (*memberIterator)->rankFitness/members.size();
                memberIterator++;
            }
            
            (*speciesIterator)->totalSharedFitness = totalFitness;
            (*speciesIterator)->speciesDist.clear();
            speciesIterator++;
        }
    }
    
    // now update the species distance maps
    for (speciesIterator = speciesList.begin();speciesIterator != speciesList.end(); speciesIterator++) {
        NEATExtendedSpecies *s1 = *speciesIterator;
        for (std::vector<NEATExtendedSpecies *>::iterator innerIter =speciesIterator+1; innerIter != speciesList.end(); innerIter++) {
            NEATExtendedSpecies *s2 = *innerIter;
            
            double d = genomeDistance(s1->representative, s2->representative);
            s1->speciesDist[s2] = d;
            s2->speciesDist[s1] = d;
        }
    }
    
}

// does s1 dominate s2?
// i.e. all s1 fitness values are greater than s2?
static int domination(SystemInfo *s1, SystemInfo *s2)
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

void MONEAT::rankSystems()
{
    std::vector<SystemInfo *>::iterator popIter;
    
    std::map<SystemInfo *, std::vector<SystemInfo *> > dominationMap;
    
    std::vector<SystemInfo *> bestFront;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        std::vector<SystemInfo *> dominated;
        SystemInfo *sys = *popIter;
        sys->dominationCount = 0;
        typename std::vector<SystemInfo *>::iterator innerIter;
        for (innerIter = population.begin(); innerIter != population.end(); innerIter++) {
            if (*innerIter == *popIter)
                continue;
            
            int dominationScore = domination( *innerIter, *popIter);
            
            if (dominationScore == -1) {
                // *popIter dominates
                dominated.push_back(*innerIter);
            }
            if (dominationScore == 1) {
                // *innerIter dominates
                sys->dominationCount += 1;
            }
        }
        if (sys->dominationCount == 0) {
            sys->rankFitness = 1.0;
            bestFront.push_back(sys);
        }
        dominationMap[sys] = dominated;
        
    }
    
    int currentRank = 1;
    std::vector<SystemInfo *> currentFront = bestFront;
    while (!currentFront.empty()) {
        std::vector<SystemInfo *> nextFront;
        std::vector<SystemInfo *>::iterator frontIterator;
        
        for (frontIterator = currentFront.begin(); frontIterator != currentFront.end(); frontIterator++) {
            std::vector<SystemInfo *> dominated = dominationMap[*frontIterator];
            std::vector<SystemInfo *>::iterator it;
            for (it = dominated.begin(); it != dominated.end(); it++) {
                SystemInfo *s = *it;
                s->dominationCount -= 1;
                // assert(s->dominationCount >=0);
                if (s->dominationCount == 0) {
                    s->rankFitness = 1.0/(currentRank+1);
                    nextFront.push_back(s);
                }
            }
        }
        
        currentRank += 1;
        currentFront = nextFront;
    }
    
    std::cout << currentRank << " ranks examined.\n";
}

// descending, not ascending, order
static bool compareSpeciesFitness( NEATExtendedSpecies *s1, const NEATExtendedSpecies *s2)
{
    return s1->totalSharedFitness > s2->totalSharedFitness;
}

bool MONEAT::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
        
        first_run = true;
    }
    
    typename std::vector<SystemInfo *>::iterator popIter;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        SystemInfo *sys = *popIter;
        sys->fitnesses.clear();
        
        typename std::vector<EvaluationFunction>::iterator funcIter;
        for (funcIter =  evaluationFunctions.begin(); funcIter != evaluationFunctions.end(); funcIter++) {
            double val = (*funcIter)(sys->individual);
            sys->fitnesses.push_back(val);
        }
    }
    
    rankSystems();
    
    bool last_run =  (generations >= maxGenerations) ;
    
#if RESTRICT_SPECIES
    {
        // clear species membership lists
        std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            (*speciesIterator)->members.clear();
            speciesIterator++;
        }
    }
#else
    speciesList.clear();
#endif
    
    speciate();
    
    updateSharedFitnesses();
    
    
    // figure out how many of each species to save
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += speciesList[i]->totalSharedFitness;
    }
    
    std::vector<SystemInfo *> individualsToSave;
    // sort the species to ensure we keep the best
    std::sort(speciesList.begin(), speciesList.end(), compareSpeciesFitness);
    
    std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo *> members = (*speciesIterator)->members;
        size_t numMembers = members.size();
        
        double  proportionToSave = ((*speciesIterator)->totalSharedFitness)/sharedFitnessSum;
        long numToSave = std::min((size_t)(proportionToSave*populationSize), numMembers);
        
        assert(numToSave >= 0);
        
        std::sort(members.begin(), members.end(), compareIndividuals);
        
        // don't let population grow unchecked
        if (individualsToSave.size() > populationSize)
            numToSave = 0;
        
        
        if (numToSave > 0) {
            std::vector<SystemInfo *> newMembers;
            std::vector<SystemInfo *>::iterator it = members.begin();
            
            while (newMembers.size() < numToSave) {
                
                SystemInfo *saved = new SystemInfo(**it);
                newMembers.push_back(saved);
                individualsToSave.push_back(saved);
                it++;
            }
            
            
            std::vector<SystemInfo *>::iterator deleteIt =  members.begin();
            
            while (deleteIt !=  members.end()) {
                delete *deleteIt;
                deleteIt =  members.erase(deleteIt);
            }
            
            (*speciesIterator)->members = newMembers;
            speciesIterator++;
            
        } else {
            speciesIterator = speciesList.erase(speciesIterator);
            
        }
        
    }
    
    
    population = individualsToSave;
    
    
    // stagnation if fitnesses are within 1% of each other
    
    
    generations++;
    
    bestIndividuals.clear();
    
    const double bestRank = 1;
    
    // find the top ranking individuals
    std::vector<SystemInfo *>::iterator it = std::find_if(individualsToSave.begin(), individualsToSave.end(), [bestRank](const SystemInfo * elem){return elem->rankFitness == bestRank;});
    for (; it != individualsToSave.end(); it = std::find_if(++it, individualsToSave.end(), [bestRank](const SystemInfo * elem){return elem->rankFitness == bestRank;})) {
        bestIndividuals.push_back(*it);
    }
    
    assert(bestIndividuals.size() > 0);
    
    
    if (!last_run) {
        spawnNextGeneration();
    }
    
    return  last_run;
}


