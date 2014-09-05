
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
#define INTERSPECIES 0

#define ROULETTE_SELECT 1
#define ELITISM 0

#define TIGHT_CLUSTERS 0
#define FINE_RANK 0


MONEAT::MONEAT(int populationSize)
:BaseNEAT(populationSize), nextSpeciesNumber(1),
stagnantGenerations(0)
{
}

void MONEAT::prepareInitialPopulation()
{
    // we need to index the genome of the 'base' system, so the genome ordering goes right
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
    
    //empty the innovations
    std::vector<InnovationInfo *>::iterator it2 = newConnections.begin();
    while (it2 != newConnections.end()) {
        delete  *it2;
        it2 = newConnections.erase(it2);
    }
    delete origin;
}

// return true if sys1 comes before sys2
static bool compareIndividuals(SystemInfo *sys1, SystemInfo *sys2)
{
    return sys1->rankFitness < sys2->rankFitness;
}
//
//static bool distanceCompare(const std::pair<NEATExtendedSpecies *, double> &p1, const std::pair<NEATExtendedSpecies*, double> &p2)
//{
//    return p1.second < p2.second;
//}

#if RESTRICT_SPECIES
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
    
    std::map<NEATExtendedSpecies *, std::vector<SystemInfo *> > newMemberMap;

    
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += (speciesList[i]->totalSharedFitness);
    }
    
    // create population children through recombination
    std::vector<SystemInfo *>newGeneration;
    
    
    std::vector<NEATExtendedSpecies *>::iterator speciesIt = speciesList.begin();
    while(speciesIt != speciesList.end()) {
        // selectionnnnn
        NEATExtendedSpecies *breedingSpecies = *speciesIt;
        std::vector<SystemInfo *>newMembers;
        
        std::vector<SystemInfo *>individuals = breedingSpecies->members;
        
#if INTERSPECIES
        std::vector<SystemInfo *> individuals2 = chooseCompatibleSpecies(breedingSpecies, maxSpeciesDist)->members;
#else
        std::vector<SystemInfo *>individuals2 = individuals;
#endif
        
        double fitness = (*speciesIt)->totalSharedFitness;
       double proportionToSave = fitness/sharedFitnessSum;

        if (speciesList.size() == 1)
            proportionToSave = 1.0;
        
        long numToSpawn = (proportionToSave*populationSize) ;//- individuals.size();
        
        if (numToSpawn < 0)
            numToSpawn = 0;
        
        while (newMembers.size() < numToSpawn) {
            
            SystemInfo *p1;
            SystemInfo *p2;
            
            MNIndividual *child;
            
#if ROULETTE_SELECT
            //stochastic acceptance
            bool accepted = false;
            while (!accepted) {
                p1 = individuals[uniformlyDistributed(individuals.size())];
                accepted = (uniformProbability() < 1.0/(p1->rankFitness));
            }
            
            
            accepted = false;
            while (!accepted) {
                p2 = individuals2[uniformlyDistributed(individuals2.size())];
                accepted = (uniformProbability() < 1.0/(p2->rankFitness));
            }
            
#else
            const int tournamentSize = 3;
  

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
        
        
        newMemberMap[*speciesIt] = newMembers;
        speciesIt++;
    }
    
    for (speciesIt = speciesList.begin(); speciesIt != speciesList.end(); speciesIt++) {
        (*speciesIt)->members = newMemberMap[*speciesIt];
    }
    
    // clear out the old pop and replace with new generation
    for (std::vector<SystemInfo *>::iterator popIter = population.begin(); popIter != population.end(); popIter++) {
        delete *popIter;
    }
    
    population = newGeneration;
    
}
#else

void MONEAT::spawnNextGeneration()
{
    // create population children through recombination
    std::vector<SystemInfo *>newGeneration;
    
    SystemInfo *p1;
    SystemInfo *p2;
    
    MNIndividual *child;
    
#if ELITISM
   // copy best individuals to this generation
    for (int i=0; i<bestIndividuals.size(); i++) {
        SystemInfo *copy = new SystemInfo(*bestIndividuals[i]);
        newGeneration.push_back(copy);
    }
#endif
        while (newGeneration.size() < populationSize) {
            
#if ROULETTE_SELECT
            //stochastic acceptance
            bool accepted = false;
            while (!accepted) {
                p1 = population[uniformlyDistributed(population.size())];
                accepted = (uniformProbability() < 1.0/(p1->rankFitness));
            }
            
            
            accepted = false;
            while (!accepted) {
                p2 = population[uniformlyDistributed(population.size())];
                accepted = (uniformProbability() < 1.0/(p2->rankFitness));
            }
            
#else
            const int tournamentSize = 3;
            
            
            // tournament
            p1 = (population[uniformlyDistributed(population.size())]);
            p2 = (population[uniformlyDistributed(population.size())]);
            
            
                for (int i=0; i<tournamentSize; i++) {
                    SystemInfo *contender = population[uniformlyDistributed(population.size())];
                    if (contender->rankFitness > p1->rankFitness)
                        p1 = contender;
                }
                
                for (int i=0; i<tournamentSize; i++) {
                    SystemInfo *contender = population[uniformlyDistributed(population.size())];
                    if (contender->rankFitness > p2->rankFitness)
                        p2 = contender;
                }
            
#endif
            child = combineSystems(p1, p2);
            // mutate the new individual (maybe)
            mutateSystem(child);
            
            SystemInfo *i = new SystemInfo(child);
            newGeneration.push_back(i);
        }
    
    // clear out the old pop and replace with new generation
    for (std::vector<SystemInfo *>::iterator popIter = population.begin(); popIter != population.end(); popIter++) {
        delete *popIter;
    }
    
    population = newGeneration;
    
}

#endif





NEATExtendedSpecies *MONEAT::chooseCompatibleSpecies(NEATExtendedSpecies *species, double maxDistance){
    NEATExtendedSpecies *chosen = speciesList[uniformlyDistributed(speciesList.size())];

    bool accepted = (maxDistance == 0); // if the distances are all 0, accept the first one we chose
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
    if (verbose)
        std::cout << nSpecies << " species (" << generations <<") \n";
    
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
    
            
#if TIGHT_CLUSTERS
            SystemInfo *centroid = members.front();
            double minDist = INFINITY;
#endif
            
            
            // update the fitnesses
            double  totalFitness = 0.0;
            std::vector<SystemInfo *>::iterator memberIterator = members.begin();
            while (memberIterator != members.end()) {
                double invFitness = 1.0/(*memberIterator)->rankFitness;
                totalFitness += invFitness/members.size();
//                double fitness = (*memberIterator)->rankFitness/members.size();
//                totalFitness += fitness;
                
#if TIGHT_CLUSTERS
                // find the central member so we can assign it as the new rep
                double myMin = INFINITY;
                for (std::vector<SystemInfo *>::iterator it = memberIterator+1; it != members.end(); it++) {
                    double d = genomeDistance((*it)->individual, (*memberIterator)->individual);
                    if (d < myMin)
                        myMin = d;
                }
                if (myMin < minDist) {
                    minDist = myMin;
                    centroid = *memberIterator;
                }
#endif
                memberIterator++;
            }
            
            (*speciesIterator)->totalSharedFitness = totalFitness;
            (*speciesIterator)->speciesDist.clear();
            
#if TIGHT_CLUSTERS
            MNIndividual *rep = centroid->individual;
#else
            long index = uniformlyDistributed(members.size());
            MNIndividual *rep = (members[index]->individual);
#endif
            
            delete (*speciesIterator)->representative;
            (*speciesIterator)->representative = rep->clone();
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
// i also consider a to dominate b if they are otherwise non dominating
// but the sum of a's fitnesses is smaller than the sum of b's
static int domination(SystemInfo *s1, SystemInfo *s2)
{
    bool firstDominates = true;
    bool secondDominates = true;
    
    std::vector<double>::iterator fitnessIter1 = s1->fitnesses.begin();
    std::vector<double>::iterator fitnessIter2 = s2->fitnesses.begin();
    
#if FINE_RANK
    double fitnessSum1 = 0.0;
    double fitnessSum2 = 0.0;
#endif
    while (fitnessIter1 != s1->fitnesses.end()) {
#if FINE_RANK
        fitnessSum1 += *fitnessIter1;
        fitnessSum2 += *fitnessIter2;
#endif
        if (*fitnessIter1 > *fitnessIter2)
            firstDominates = false;
        
        if (*fitnessIter2 > *fitnessIter1)
            secondDominates = false;
        fitnessIter2++;
        fitnessIter1++;
    }
    
    if (firstDominates == secondDominates) {
#if FINE_RANK
        if (fitnessSum1 < fitnessSum2)
            return 1;
        if (fitnessSum2 < fitnessSum1)
            return -1;
#endif
        return 0;
    }
    
    if (firstDominates)
        return 1;
    
    if (secondDominates)
        return -1;
    
    return 0;
}

// using NSGA-II
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
                    s->rankFitness = (currentRank+1);
                    nextFront.push_back(s);
                }
            }
        }
        
        currentRank += 1;
        currentFront = nextFront;
    }
    if (verbose)
        std::cout << currentRank << " ranks examined.\n";
}


void MONEAT::tick()
{
    if (population.size() == 0) {
        prepareInitialPopulation();
    } else {
        spawnNextGeneration();
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
    
    
    
#if RESTRICT_SPECIES
    bool respeciate = first_run;
    if (!respeciate) {
        respeciate = (generations % 25 == 0);
    }
#else
    bool respeciate = true;
#endif

    if (respeciate) {
        
        // clear species membership lists
        std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            (*speciesIterator)->members.clear();
            speciesIterator++;
        }
        
        speciate();
    }

    
    updateSharedFitnesses();
    
    std::sort(population.begin(), population.end(), compareIndividuals);

    std::map<NEATExtendedSpecies *, long> preservationMap;
    
    // figure out how many of each species to save
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += (speciesList[i]->totalSharedFitness);
    }
    
    std::vector<NEATExtendedSpecies *>::iterator speciesIterator = speciesList.begin();
    
    while(speciesIterator != speciesList.end()) {
        
        // for debug
       // std::sort((*speciesIterator)->members.begin(), (*speciesIterator)->members.end(), compareIndividuals);
        
        double fitness = (*speciesIterator)->totalSharedFitness;        
        double proportionToSave = (fitness)/sharedFitnessSum;
        long numToSave = std::min((size_t)(proportionToSave*populationSize), (*speciesIterator)->members.size());
        
//        if (numToSave > 1)
//            numToSave /= 2;
        
        if (numToSave > 0) {
            std::sort((*speciesIterator)->members.begin(), (*speciesIterator)->members.end(), compareIndividuals);
            preservationMap[(*speciesIterator)] = numToSave;
            speciesIterator++;
        } else {
            delete (*speciesIterator);
            speciesIterator = speciesList.erase(speciesIterator);
        }
    }
    
    std::vector<SystemInfo *>individualsToSave;
    std::map<NEATExtendedSpecies *, std::vector<SystemInfo *> > newMemberMap;
    
    bool done;
    long i=0;
    do {
        done = true;
        std::vector<NEATExtendedSpecies *>::iterator speciesIterator;
        for (speciesIterator = speciesList.begin(); speciesIterator != speciesList.end(); speciesIterator++) {
            if (preservationMap[*speciesIterator] > i) {
                std::vector<SystemInfo *> members = (*speciesIterator)->members;
                SystemInfo *saved = new SystemInfo(*(*(members.begin() + i)));
                individualsToSave.push_back(saved);

                newMemberMap[*speciesIterator].push_back(saved);
                done = false;
            }
        }
        i++;
    } while(!done);
    
    // replace the species member lists
    for (speciesIterator = speciesList.begin(); speciesIterator != speciesList.end(); speciesIterator++) {
        (*speciesIterator)->members = newMemberMap[*speciesIterator];
    }
    
    std::sort(individualsToSave.begin(), individualsToSave.end(), compareIndividuals);

    
    bestIndividuals.clear();
    
    const double bestRank = 1;
    
    // find the top ranking individuals
   std::vector<SystemInfo *>::iterator it = std::find_if(individualsToSave.begin(), individualsToSave.end(), [bestRank](const SystemInfo * elem){return elem->rankFitness == bestRank;});
    for (; it != individualsToSave.end(); it = std::find_if(++it, individualsToSave.end(), [bestRank](const SystemInfo * elem){return elem->rankFitness == bestRank;})) {
        bestIndividuals.push_back(*it);
    }
    
    assert(bestIndividuals.size() > 0);
    
    // delete the old population - saved individuals have been copied
    it = population.begin();
    while (it != population.end()) {
        delete *it;
        it = population.erase(it);
    }
    
    population = individualsToSave;

    generations++;

    
if (verbose)
    std::cout << "\n";
    
}


