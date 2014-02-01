
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

#define RESTRICT_SPECIES 0
#define ROULETTE_SELECT 0
#define INTERSPECIES 0

template <class IndividualType, class InnovationType>
MONEAT <IndividualType, InnovationType>::MONEAT(int populationSize, int maxGenerations, int maxStagnation)
:populationSize(populationSize), maxStagnation(maxStagnation), maxGenerations(maxGenerations)

{
    _bestIndividual = NULL;
    nextInnovationNumber = 1;
    nextSpeciesNumber = 1;
    allTimeBestFitness = -INFINITY;
    generations = 0;
    stagnantGenerations = 0;
    origin = NULL;
}

template <class IndividualType, class InnovationType>
MONEAT <IndividualType, InnovationType>::~MONEAT()
{
    //empty the species lists
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo<IndividualType> *> members = speciesIterator->members;
        typename std::vector<SystemInfo<IndividualType> *>::iterator it = members.begin();
        while (it != members.end()) {
            it = members.erase(it);
        }
        speciesIterator = speciesList.erase(speciesIterator);
    }
    delete origin;
}

// return true if sys1 comes before sys2
template <class IndividualType>
static bool compareIndividuals(SystemInfo<IndividualType> *sys1, SystemInfo<IndividualType> *sys2)
{
    return sys1->rankFitness > sys2->rankFitness;
}

template <class InnovationType>
static bool compareInnovationNumbers(const InnovationType a1, const InnovationType a2)
{
    return a1.innovationNumber < a2.innovationNumber;
}

template <class IndividualType, class InnovationType>
IndividualType *MONEAT <IndividualType, InnovationType>::combineSystems(SystemInfo<IndividualType> *sys1, SystemInfo<IndividualType> *sys2)
{
    
    if (sys1 == sys2) {
        // breeding with myself?
        return new IndividualType(*sys1->individual);
    }
    
    IndividualType *newChild = createInitialIndividual();
    std::vector<InnovationType> genome1 = sys1->individual->connectionGenome();
    std::vector<InnovationType> genome2 = sys2->individual->connectionGenome();
    
    // holy crap, c++ has some good standard methods, like *set difference*
    std::vector<InnovationType> matchingGenes1;
    std::vector<InnovationType> matchingGenes2;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers<InnovationType>);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers<InnovationType>);
    
    
    
    // first we handle the matching genes
    for (int i=0; i<matchingGenes1.size(); i++){
        InnovationType gene = matchingGenes1[i];
        IndividualType *parent = sys1->individual;
        // randomly choose 1 to add to the child
        double selector = (double )rand()/RAND_MAX;
        
        if (selector <= p_c) {
            parent = sys2->individual;
            gene = matchingGenes2[i];
        }
        newChild->addGeneFromParentSystem(*parent, gene);
    }
    
    
    std::vector<InnovationType> disjointandExcess1;
    std::vector<InnovationType> disjointandExcess2;
    
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointandExcess1), compareInnovationNumbers<InnovationType>);
    std::set_difference(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(), std::back_inserter(disjointandExcess2), compareInnovationNumbers<InnovationType>);
    
    
    // then the disjoint genes
    if (sys1->rankFitness >= sys2->rankFitness) {
        for (int i=0; i<disjointandExcess1.size(); i++) {
            newChild->addGeneFromParentSystem(*sys1->individual, disjointandExcess1[i]);
        }
    }
    
    // if rankFitness is equal, include all genes
    if (sys2->rankFitness >= sys1->rankFitness) {
        for (int i=0; i<disjointandExcess2.size(); i++) {
            newChild->addGeneFromParentSystem(*sys2->individual, disjointandExcess2[i]);
        }
    }
    return newChild;
}

template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::spawnNextGeneration()
{
    double  rankFitnessSum = 0.0;
    typename std::vector<SystemInfo<IndividualType> *>::iterator it = population.begin();
    while (it != population.end()) {
        rankFitnessSum += (*it)->rankFitness;
        it++;
    }

    // create population children through recombination
    std::vector<SystemInfo<IndividualType> *>newGeneration;
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIt = speciesList.begin();
    while(speciesIt != speciesList.end()) {
        // selectionnnnn
        NEATSpecies<IndividualType> breedingSpecies = *speciesIt;//chooseBreedingSpecies(rankFitnessSum);
        std::vector<SystemInfo<IndividualType> *>newMembers;
        
        std::vector<SystemInfo<IndividualType> *>individuals = breedingSpecies.members;
        
        while (newMembers.size() < individuals.size()) {
            SystemInfo<IndividualType> *p1;
            SystemInfo<IndividualType> *p2;
            
            IndividualType *child;
            
            int parentPosition = arc4random_uniform((int)individuals.size());
            p1 = (individuals[parentPosition]);
#if INTERSPECIES
            double interspecies_mate = (double)rand()/RAND_MAX;
            if (stagnantGenerations > maxStagnation/3 || interspecies_mate < 0.1) {
                breedingSpecies = chooseBreedingSpecies(rankFitnessSum);
                individuals = breedingSpecies.members;
                
                int parentPosition = arc4random_uniform((int)individuals.size());
                p2 = (individuals[parentPosition]);
            } else {
#endif
                if (individuals.size() > 1) {
                    
                    int parentPosition2;
                    do {
                        parentPosition2 = arc4random_uniform((int)individuals.size());
                    } while (parentPosition2 == parentPosition);
                    p2 = (individuals[parentPosition2]);
                } else {
                    p2 = p1;
                }
#if INTERSPECIES
            }
#endif
            
            child = combineSystems(p1, p2);
            // mutate the new individual (maybe)
            mutateSystem(child);
            
            SystemInfo<IndividualType> *i = new SystemInfo<IndividualType>(child);
            newGeneration.push_back(i);
            newMembers.push_back(i);
        }
        speciesIt->members.insert(speciesIt->members.end(), newMembers.begin(), newMembers.end());
        speciesIt++;
    }
    // add the original population too
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());

//    population = newGeneration;
}

template <class IndividualType, class InnovationType>
NEATSpecies<IndividualType> &MONEAT<IndividualType, InnovationType>::chooseBreedingSpecies(double totalFitness)
{
    
    // choose a species according to species fitness
    // then chose parents at random from within that species
    
    // Assume sorted population
    double selector = ((double )rand()/RAND_MAX);
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator it = speciesList.begin();
    double  cumulativeFitness = 0.0;
    while (it != speciesList.end()) {
        cumulativeFitness += it->totalSharedFitness/totalFitness;
        if (cumulativeFitness >= selector) {
            break;
        }
        it++;
    }
    
    if (it == speciesList.end())
        it--;
    
    return *it;
}


template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::speciate()
{
    typename std::vector<SystemInfo<IndividualType> *>::iterator populationIter = population.begin();
    
    while (populationIter != population.end()) {
        // assign to species
        bool added= false;
        
        typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            double  distance = genomeDistance((*populationIter)->individual, speciesIterator->representative);
            if (fabs(distance) < d_threshold) {
                added = true;
                speciesIterator->members.push_back((*populationIter));
                break;
            }
            speciesIterator++;
        }
        
        if (!added) {
            // new species!!
            NEATSpecies<IndividualType> s = NEATSpecies<IndividualType>();
            s.representative = new IndividualType(*(*populationIter)->individual);
            s.members.push_back(*populationIter);
            speciesList.push_back(s);
            s.speciesNumber = nextSpeciesNumber++;
        }
        
        populationIter++;
    }
    
    
}

template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::updateSharedFitnesses()
{
    long nSpecies = speciesList.size();
    fprintf(stderr, "%ld species\n", nSpecies);
    
    // for each species, copy a random current member to be the representative for the next generation, and adjust the fitnesses for sharing
    // kill off a species with no members
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        // pick a new rep
        std::vector<SystemInfo<IndividualType> *> members = speciesIterator->members;
        if (members.size() == 0) {
            NEATSpecies<IndividualType> extinctSpecies = *speciesIterator;
            speciesIterator = speciesList.erase(speciesIterator);
        } else {
            long index = uniformlyDistributed(members.size());
            IndividualType *rep = (members[index]->individual);
            delete speciesIterator->representative;
            speciesIterator->representative = new IndividualType(*rep);
            
            // fprintf(stderr, "\tSpecies %d: %ld members\n", ++i, members.size());
            
            // update the fitnesses
            double  totalFitness = 0.0;
            typename std::vector<SystemInfo<IndividualType> *>::iterator memberIterator = members.begin();
            while (memberIterator != members.end()) {
                assert((*memberIterator)->rankFitness != -INFINITY);
                (*memberIterator)->rankFitness /= members.size();
                totalFitness += (*memberIterator)->rankFitness;
                memberIterator++;
            }
            
            speciesIterator->totalSharedFitness = totalFitness;
            speciesIterator++;
        }
    }
}

// unlike base NEAT, also take node differences into account
template <class IndividualType, class InnovationType>
double  MONEAT <IndividualType, InnovationType>::genomeDistance( IndividualType *sys1,  IndividualType *sys2)
{
    std::vector<InnovationType> genome1 = sys1->connectionGenome();
    std::vector<InnovationType> genome2 = sys2->connectionGenome();
    
    // find the longer one
    size_t longerSize = (genome1.size() > genome2.size() ? genome1.size()  : genome2.size());
    
    size_t nDisjoint = 0;
    size_t nExcess = 0;
    size_t nMatching = 0;
    size_t matchingDiff = 0;
    
    // yay stdlib!
    std::vector<InnovationType> matchingGenes1;
    std::vector<InnovationType> matchingGenes2;
    std::vector<InnovationType> disjointFromSys1;
    std::vector<InnovationType> disjointFromSys2;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers<InnovationType>);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers<InnovationType>);
    
    // difference takes things in the first range that are not in the second - we will have to check for excess ourself
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointFromSys1), compareInnovationNumbers<InnovationType>);
    std::set_difference(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(), std::back_inserter(disjointFromSys2), compareInnovationNumbers<InnovationType>);
    
    
    nMatching = matchingGenes2.size();
    // first find the distance between matching genes
    for (int i=0; i<nMatching; i++) {
        InnovationType a1 = matchingGenes1[i];
        InnovationType a2 = matchingGenes2[i];
        
        matchingDiff += IndividualType::connectionDifference(a1, a2);
    }
    
    // now determine the excess vs disjoint
    // if one of the disjoint sets is empty, then it's all excess
    if (disjointFromSys2.size() == 0 || disjointFromSys1.size() == 0) {
        nExcess = disjointFromSys1.size() + disjointFromSys2.size();
    } else {
        // else chalk everything up to disjoint
        nDisjoint = disjointFromSys1.size() + disjointFromSys2.size();
        
        // then find the number of excess genes
        int last1 = disjointFromSys1.back().innovationNumber;
        int last2 = disjointFromSys2.back().innovationNumber;
        
        std::vector<InnovationType> shorter;
        std::vector<InnovationType> longer;
        
        if (last1 != last2) {
            if (last1 > last2) {
                // genome 1 has excess
                shorter = disjointFromSys2;
                longer = disjointFromSys1;
            } else if (last1 < last2) {
                // genome 2 has excess, switch last1 and last2
                shorter=disjointFromSys1;
                longer = disjointFromSys2;
                
                int temp = last1;
                last1 = last2;
                last2 = temp;
            }
            // go back from end of longer until we run into a value <= the end of shorter
            size_t i;
            for (i=longer.size(); i; i--) {
                if (longer[i].innovationNumber < shorter.back().innovationNumber)
                    break;
            }
            nExcess = i;
            nDisjoint -= nExcess;
        }
    }
    
    
    double diff_node = sys1->nodeDifference(*sys2);
    long moreNodes = std::min(sys1->numberOfNodes(), sys2->numberOfNodes());
    
    double  d = w_disjoint*nDisjoint/longerSize + w_excess*nExcess/longerSize + w_matching*matchingDiff/nMatching;
    
    //    double  d = w_disjoint*nDisjoint + w_excess*nExcess + w_matching*matchingDiff;
    
    //  d+= + w_matching_node*diff_node/moreNodes;
    
    //assert(!isUnreasonable(d));
    return d;
}


template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::prepareInitialPopulation()
{
    // mutate the initial system to get an initial population
    while (population.size() < populationSize) {
        IndividualType *newSystem = createInitialIndividual();
        mutateSystem(newSystem);
        SystemInfo<IndividualType> *info = new SystemInfo<IndividualType>(newSystem);
        population.push_back(info);
    }
}

// does s1 dominate s2?
// i.e. all s1 fitness values are greater than s2?
template <class IndividualType>
static int domination(SystemInfo<IndividualType> *s1, SystemInfo<IndividualType> *s2)
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
    
    if (firstDominates)
        return 1;
    else if (secondDominates)
        return -1;
    
    return 0;
}

template <class IndividualType, class InnovationType>
std::vector<SystemInfo<IndividualType> *> MONEAT <IndividualType, InnovationType>::rankSystems()
{
    typename std::vector<SystemInfo<IndividualType> *>::iterator popIter;
    
    std::map<SystemInfo<IndividualType> *, std::vector<SystemInfo<IndividualType> *> > dominationMap;
    
    
    std::vector<SystemInfo<IndividualType> *> bestFront;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        std::vector<SystemInfo<IndividualType> *> dominated;
        SystemInfo<IndividualType> *sys = *popIter;
        sys->dominationCount = 0;
        
        typename std::vector<SystemInfo<IndividualType> *>::iterator innerIter;
        for (innerIter = population.begin(); innerIter != population.end(); innerIter++) {
            if (*innerIter == *popIter)
                continue;
            
            int dominationScore = domination(*popIter, *innerIter);
            
            if (dominationScore == 1) {
                // *popIter dominates
                dominated.push_back(*innerIter);
            }
            if (dominationScore == -1) {
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
    
    std::vector<SystemInfo<IndividualType> *> currentFront = bestFront;
    
    while (!currentFront.empty()) {
        std::vector<SystemInfo<IndividualType> *> nextFront;
        
        typename std::vector<SystemInfo<IndividualType> *>::iterator frontIterator;
        for (frontIterator = currentFront.begin(); frontIterator != currentFront.end(); frontIterator++) {
            typename std::vector<SystemInfo<IndividualType> *>::iterator it;
            for (it = dominationMap[*frontIterator].begin(); it != dominationMap[*frontIterator].end(); it++) {
                SystemInfo<IndividualType> *s = *it;
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
    return bestFront;
}

// descending, not ascending, order
template <class IndividualType>
static bool compareSpeciesFitness(const NEATSpecies<IndividualType>& s1, const NEATSpecies<IndividualType>& s2)
{
    return s1.totalSharedFitness > s2.totalSharedFitness;
}

template <class IndividualType, class InnovationType>
bool MONEAT <IndividualType, InnovationType>::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
        
        origin = new IndividualType(*population.front()->individual);
        first_run = true;
    }
    
    
    typename std::vector<EvaluationFunction>::iterator funcIter;
    for (funcIter =  evaluationFunctions.begin(); funcIter != evaluationFunctions.end(); funcIter++) {
        typename std::vector<SystemInfo<IndividualType> *>::iterator popIter;
        for (popIter = population.begin(); popIter != population.end(); popIter++) {
            SystemInfo<IndividualType> *sys = *popIter;
            double val = (*funcIter)(sys->individual);
            sys->fitnesses.push_back(val);
        }
    }
    
    std::vector<SystemInfo<IndividualType> *> bestSystems = rankSystems();
    
    bool last_run =  (generations >= maxGenerations) ;
    if (last_run) {
        SystemInfo<IndividualType> *someone = bestSystems.front();
        std::cout << bestSystems.size() << " systems in best rank.\n";
    }
    
#if !RESTRICT_SPECIES
   {
       // clear species membership lists
       typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
       while (speciesIterator != speciesList.end()) {
           speciesIterator->members.clear();
           speciesIterator++;
       }
   }
#endif
    
#if RESTRICT_SPECIES
    if (first_run) {
#endif
        speciate();
#if RESTRICT_SPECIES
    } else {
    
          bool regroup = false;
          if (speciesList.size() < 10 || stagnantGenerations > maxStagnation/5) {
              regroup = true;
              d_threshold /= 1.05;
          }
      
          if (speciesList.size() > 20) {
              regroup = true;
              d_threshold *= 1.05;
          }
      
          if (regroup) {
              // empty the species member lists and remake them
              speciesList.clear();
              speciate();
          }
          
    }
#endif
    updateSharedFitnesses();

    
    double  bestFitness = -INFINITY; // we want zero fitness
    // figure out how many of each species to save
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += speciesList[i].totalSharedFitness;
    }
    
    std::vector<SystemInfo<IndividualType> *> individualsToSave;
    // sort the species to ensure we keep the best
    std::sort(speciesList.begin(), speciesList.end(), compareSpeciesFitness<IndividualType>);
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    IndividualType *bestFound = NULL;
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo<IndividualType> *> members =speciesIterator->members;
        size_t numMembers = members.size();

        double  proportionToSave = (speciesIterator->totalSharedFitness)/sharedFitnessSum;
        long numToSave = std::min((size_t)(proportionToSave*populationSize), numMembers);
        
        assert(numToSave >= 0);
        
        std::sort(members.begin(), members.end(), compareIndividuals<IndividualType>);
        
        // don't let population grow unchecked
        if (individualsToSave.size() > populationSize)
            numToSave = 0;
        
        
        //stochastic universal selection
        if (numToSave > 0) {
            std::vector<SystemInfo<IndividualType> *> newMembers;
            typename std::vector<SystemInfo<IndividualType> *>::iterator it = members.begin();

#if ROULETTE_SELECT
            double choiceSep = (speciesIterator->totalSharedFitness/numToSave);
            double pointer = ((double)rand()/RAND_MAX)*choiceSep;
            std::vector<long> savedPositions;
            double cumFitness = (*it)->rankFitness;
            while (newMembers.size() < numToSave) {
                if (cumFitness >= pointer) {
                    pointer += choiceSep;
                    
                    SystemInfo<IndividualType> *individual = *it;

                    IndividualType *copied = new IndividualType(*individual->individual);
                    
                    SystemInfo<IndividualType> *saved = new SystemInfo<IndividualType>(copied);
                    double  rankFitness = individual->rankFitness * numMembers;
                    saved->rankFitness = rankFitness;
                    newMembers.push_back(saved);
                    individualsToSave.push_back(saved);
                    if (rankFitness > bestFitness) {
                        bestFitness = rankFitness;
                        bestFound = copied;
                    }
                    
                } else {
                    it++;
                    assert(it != members.end());
                    cumFitness += (*it)->rankFitness;
                }
            }
#else
            while (newMembers.size() < numToSave) {
                
                SystemInfo<IndividualType> *individual = *it;
                
                IndividualType *copied = new IndividualType(*individual->individual);
                
                SystemInfo<IndividualType> *saved = new SystemInfo<IndividualType>(copied);
                double  rankFitness = individual->rankFitness * numMembers;
                saved->rankFitness = individual->rankFitness;
                newMembers.push_back(saved);
                individualsToSave.push_back(saved);
                if (rankFitness > bestFitness) {
                    bestFitness = rankFitness;
                    bestFound = copied;
                }
                it++;
            }
#endif
            
            typename std::vector<SystemInfo<IndividualType> *>::iterator deleteIt =  members.begin();
            
            while (deleteIt !=  members.end()) {
                assert(*deleteIt != NULL);
                delete *deleteIt;
                deleteIt =  members.erase(deleteIt);
            }
            
            speciesIterator->members = newMembers;
            speciesIterator++;

        } else {
            speciesIterator = speciesList.erase(speciesIterator);
 
        }
        
    }
    if (_bestIndividual)
        delete _bestIndividual;
    
    _bestIndividual = new IndividualType(*bestFound);
    
    population = individualsToSave;
    
    if (bestFitness > allTimeBestFitness)
        allTimeBestFitness = bestFitness;
    
    // stagnation if fitnesses are within 1% of each other
    
    if (fabs( 1 - (lastBestFitness/bestFitness)) < .01)
        stagnantGenerations++;
    else
        stagnantGenerations = 0;
    
    generations++;
    
    lastBestFitness = bestFitness;
   
    if (!last_run) {
        spawnNextGeneration();
    }
    
    return  last_run;
}

template <class IndividualType, class InnovationType>
IndividualType  *MONEAT <IndividualType, InnovationType>::bestIndividual()
{
    return _bestIndividual;
}

template <class IndividualType, class InnovationType>
long MONEAT <IndividualType, InnovationType>::getNumberOfIterations()
{
    return generations;
}


template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::mutateSystem(IndividualType *original)
{
    
    std::vector<InnovationType> allAttachments = original->connectionGenome();
    typename std::vector<InnovationType>::iterator it = allAttachments.begin();
    while (it != allAttachments.end()) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_conn)
            original->mutateConnectionWeight();
        it++;
    }
    
    long nNodes = original->numberOfNodes();
    for (long i=0; i<nNodes; i++) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_node)
            original->mutateNode(i);
    }
    
    // add attachment or add node - or both
    double selector1 = (double )rand()/RAND_MAX;
    if (selector1 < p_m_conn_ins) {
        // insert a new connection
        // update the innovation number
        InnovationType created = original->createConnection();
        assignInnovationNumberToAttachment(original, created);
    }
    
    double selector2 = (double )rand()/RAND_MAX;
    if (selector2 < p_m_node_ins) {
        // add a node somewhere if possible
        std::vector<InnovationType> new_edges = original->insertNode();
        for (int i=0; i<new_edges.size(); i++) {
            assignInnovationNumberToAttachment(original,  new_edges[i]);
        }
    }
    
}

template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::assignInnovationNumberToAttachment(IndividualType *individual, InnovationType i){
    // have we already created this 'innovation' in this generation?
    typename std::vector<InnovationType>::iterator it = std::find(newConnections.begin(), newConnections.end(), i);
    if (it!= newConnections.end()) {
        i.innovationNumber = it->innovationNumber;
    } else {
        i.innovationNumber = nextInnovationNumber++;
        newConnections.push_back(i);
    }
    individual->updateInnovationNumber(i);
}

template <class IndividualType, class InnovationType>
void MONEAT <IndividualType, InnovationType>::logPopulationStatistics()
{
   
}
